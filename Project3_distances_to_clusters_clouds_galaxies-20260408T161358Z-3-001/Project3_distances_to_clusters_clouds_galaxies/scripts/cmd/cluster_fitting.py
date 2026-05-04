from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

# -----------------------------
# Load data
# -----------------------------
pleiades = pd.read_csv(ROOT / "gaia" / "gaia_pleiades.csv")
hyades_ms = pd.read_csv(ROOT / "gaia" / "hyades_members_absolute_cmd.csv")

# Clean
pleiades = pleiades.dropna(subset=["ra", "dec", "parallax", "pmra", "pmdec",
                                   "phot_g_mean_mag", "bp_rp", "ruwe"]).copy()
pleiades = pleiades[(pleiades["parallax"] > 0) & (pleiades["ruwe"] < 1.4)].copy()

hyades_ms = hyades_ms.dropna(subset=["bp_rp", "M_G"]).copy()

# -----------------------------
# 1. Pleiades VPD
# Pleiades approx:
# distance ~ 135 pc => parallax ~ 7.4 mas
# pmra ~ +19 mas/yr, pmdec ~ -45 mas/yr
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(pleiades["pmra"], pleiades["pmdec"], s=5, alpha=0.35)
plt.xlabel(r"$\mu_{\alpha *}$ = pmRA [mas/yr]")
plt.ylabel(r"$\mu_\delta$ = pmDec [mas/yr]")
plt.title("Pleiades: Vector Point Diagram")
plt.grid(True)
plt.savefig(PLOT_DIR / "pleiades_vpd_all.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 2. Select likely Pleiades members
# Broad but clean first pass
# -----------------------------
members = pleiades[
    (pleiades["pmra"] > 10) & (pleiades["pmra"] < 30) &
    (pleiades["pmdec"] > -60) & (pleiades["pmdec"] < -30) &
    (pleiades["parallax"] > 6.0) & (pleiades["parallax"] < 9.0)
].copy()

field = pleiades.loc[~pleiades.index.isin(members.index)].copy()

members["distance_pc"] = 1000 / members["parallax"]

print("Total Pleiades field stars:", len(pleiades))
print("Selected Pleiades members:", len(members))
print("\nPleiades parallax distance statistics:")
print(members["distance_pc"].describe())
print("\nMedian Pleiades parallax distance:")
print(np.median(members["distance_pc"]), "pc")

members.to_csv(ROOT / "gaia" / "pleiades_members_selected.csv", index=False)

# -----------------------------
# 3. VPD selected members
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(field["pmra"], field["pmdec"], s=4, alpha=0.18, label="Field")
plt.scatter(members["pmra"], members["pmdec"], s=12, alpha=0.9, label="Selected Pleiades")
plt.xlabel(r"$\mu_{\alpha *}$ = pmRA [mas/yr]")
plt.ylabel(r"$\mu_\delta$ = pmDec [mas/yr]")
plt.title("Pleiades: VPD Member Selection")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "pleiades_vpd_members.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 4. Pleiades apparent CMD
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(field["bp_rp"], field["phot_g_mean_mag"], s=4, alpha=0.15, label="Field")
plt.scatter(members["bp_rp"], members["phot_g_mean_mag"], s=12, alpha=0.9, label="Selected Pleiades")
plt.gca().invert_yaxis()
plt.xlabel("BP - RP")
plt.ylabel("G")
plt.title("Pleiades CMD after VPD + Parallax Selection")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "pleiades_cmd_members.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 5. Prepare Hyades calibrated main sequence
# Use color range where both sequences overlap
# -----------------------------
# Sort Hyades by color
hyades_ms = hyades_ms.sort_values("bp_rp")

# Remove extreme outliers for a cleaner reference
hyades_ref = hyades_ms[
    (hyades_ms["bp_rp"] > 0.4) &
    (hyades_ms["bp_rp"] < 2.8) &
    (hyades_ms["M_G"] > 2) &
    (hyades_ms["M_G"] < 13)
].copy()

pleiades_fit = members[
    (members["bp_rp"] > hyades_ref["bp_rp"].min()) &
    (members["bp_rp"] < hyades_ref["bp_rp"].max())
].copy()

# Bin Hyades MS to make a smoother reference curve
bins = np.linspace(hyades_ref["bp_rp"].min(), hyades_ref["bp_rp"].max(), 35)
hyades_ref["bin"] = pd.cut(hyades_ref["bp_rp"], bins=bins)

hyades_binned = hyades_ref.groupby("bin", observed=True).agg(
    bp_rp=("bp_rp", "median"),
    M_G=("M_G", "median")
).dropna()

# Interpolation function: color -> absolute magnitude
f_ms = interp1d(
    hyades_binned["bp_rp"],
    hyades_binned["M_G"],
    bounds_error=False,
    fill_value=np.nan
)

# -----------------------------
# 6. Fit distance modulus
# We shift Hyades M_G by mu and match Pleiades G:
# G_pleiades ≈ M_G_hyades + mu
# -----------------------------
def chi2(mu):
    pred_G = f_ms(pleiades_fit["bp_rp"].values) + mu
    obs_G = pleiades_fit["phot_g_mean_mag"].values
    mask = np.isfinite(pred_G)
    return np.nanmedian((obs_G[mask] - pred_G[mask])**2)

result = minimize_scalar(chi2, bounds=(4.0, 8.0), method="bounded")
mu_best = result.x

d_ms_pc = 10 ** ((mu_best + 5) / 5)

print("\nBest-fit distance modulus from MS fitting:")
print(mu_best)

print("\nPleiades distance from MS fitting:")
print(d_ms_pc, "pc")

print("\nPleiades median Gaia parallax distance:")
print(np.median(members["distance_pc"]), "pc")

# Save result
with open(ROOT / "gaia" / "pleiades_ms_fit_result.txt", "w") as f:
    f.write(f"Best-fit distance modulus mu = {mu_best:.4f}\n")
    f.write(f"MS-fitting distance = {d_ms_pc:.2f} pc\n")
    f.write(f"Median Gaia parallax distance = {np.median(members['distance_pc']):.2f} pc\n")
    f.write(f"Number of selected Pleiades members = {len(members)}\n")

# -----------------------------
# 7. Overlay Hyades calibrated MS shifted onto Pleiades CMD
# -----------------------------
plt.figure(figsize=(7, 6))

plt.scatter(members["bp_rp"], members["phot_g_mean_mag"],
            s=14, alpha=0.75, label="Pleiades members")

plt.plot(hyades_binned["bp_rp"], hyades_binned["M_G"] + mu_best,
         linewidth=2.5, label=f"Hyades MS shifted by μ = {mu_best:.2f}")

plt.gca().invert_yaxis()
plt.xlabel("BP - RP")
plt.ylabel("G")
plt.title("Main Sequence Fitting: Hyades MS → Pleiades")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "pleiades_ms_fit_hyades_overlay.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 8. Show absolute Hyades MS and apparent Pleiades together
# -----------------------------
plt.figure(figsize=(7, 6))

plt.scatter(hyades_ref["bp_rp"], hyades_ref["M_G"],
            s=10, alpha=0.45, label="Hyades absolute MS")

plt.scatter(members["bp_rp"], members["phot_g_mean_mag"],
            s=10, alpha=0.65, label="Pleiades apparent CMD")

plt.gca().invert_yaxis()
plt.xlabel("BP - RP")
plt.ylabel("Magnitude")
plt.title("Hyades Absolute MS vs Pleiades Apparent CMD")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "hyades_absolute_vs_pleiades_apparent.png", dpi=300, bbox_inches="tight")
plt.show()

print("\nSaved plots in:", PLOT_DIR)