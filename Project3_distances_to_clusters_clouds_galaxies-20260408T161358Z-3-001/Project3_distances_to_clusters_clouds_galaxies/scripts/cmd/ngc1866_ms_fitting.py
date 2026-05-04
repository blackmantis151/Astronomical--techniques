from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

# -----------------------------
# Load NGC 1866 BV photometry
# -----------------------------
ngc_file = ROOT / "ngc1866_ms" / "J_AJ_118_2839_stars.csv"
ngc = pd.read_csv(ngc_file)

print("NGC 1866 columns:")
print(ngc.columns)
print(ngc.head())

ngc = ngc.dropna(subset=["Vmag", "B-V"]).copy()
ngc = ngc[(ngc["Vmag"] > 10) & (ngc["Vmag"] < 24)]
ngc = ngc[(ngc["B-V"] > -0.5) & (ngc["B-V"] < 2.0)]

# -----------------------------
# Plot raw NGC 1866 CMD
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(ngc["B-V"], ngc["Vmag"], s=2, alpha=0.35)
plt.gca().invert_yaxis()
plt.xlabel("B - V")
plt.ylabel("V")
plt.title("NGC 1866 CMD")
plt.grid(True)
plt.savefig(PLOT_DIR / "ngc1866_cmd_raw.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Approximate reference MS from Hyades Gaia CMD
# WARNING:
# This is approximate because Gaia BP-RP/M_G is not the same as Johnson B-V/M_V.
# For a course-level consistency check, we use:
# B-V ~ 0.75*(BP-RP)
# M_V ~ M_G
# -----------------------------
hyades = pd.read_csv(ROOT / "gaia" / "hyades_members_absolute_cmd.csv")
hyades = hyades.dropna(subset=["bp_rp", "M_G"]).copy()

hyades["B_V_approx"] = 0.75 * hyades["bp_rp"]
hyades["M_V_approx"] = hyades["M_G"]

hyades_ref = hyades[
    (hyades["B_V_approx"] > -0.1) &
    (hyades["B_V_approx"] < 1.5) &
    (hyades["M_V_approx"] > 1.0) &
    (hyades["M_V_approx"] < 9.0)
].copy()

hyades_ref = hyades_ref.sort_values("B_V_approx")

# Bin Hyades reference sequence
bins = np.linspace(hyades_ref["B_V_approx"].min(), hyades_ref["B_V_approx"].max(), 30)
hyades_ref["bin"] = pd.cut(hyades_ref["B_V_approx"], bins=bins)

hy_binned = hyades_ref.groupby("bin", observed=True).agg(
    B_V=("B_V_approx", "median"),
    M_V=("M_V_approx", "median")
).dropna()

f_ms = interp1d(
    hy_binned["B_V"],
    hy_binned["M_V"],
    bounds_error=False,
    fill_value=np.nan
)

# -----------------------------
# Select a plausible NGC 1866 MS region
# The cluster is young, so the bright upper MS is useful.
# Avoid red giants and field contamination.
# -----------------------------
fit_data = ngc[
    (ngc["B-V"] > -0.1) &
    (ngc["B-V"] < 0.8) &
    (ngc["Vmag"] > 16.0) &
    (ngc["Vmag"] < 23.0)
].copy()

# -----------------------------
# Fit distance modulus + rough color offset
# NGC 1866 is reddened and uses different filters.
# We fit only vertical shift first; then test color offset grid.
# -----------------------------
def chi2_for(mu, color_shift):
    color_corr = fit_data["B-V"].values - color_shift
    pred_V = f_ms(color_corr) + mu
    obs_V = fit_data["Vmag"].values
    mask = np.isfinite(pred_V)
    if mask.sum() < 20:
        return np.inf
    return np.nanmedian((obs_V[mask] - pred_V[mask])**2)

best = {"chi2": np.inf, "mu": None, "color_shift": None}

# color_shift roughly E(B-V); try 0.00 to 0.30
for color_shift in np.linspace(0.0, 0.30, 61):
    res = minimize_scalar(
        lambda mu: chi2_for(mu, color_shift),
        bounds=(17.5, 19.5),
        method="bounded"
    )
    if res.fun < best["chi2"]:
        best = {"chi2": res.fun, "mu": res.x, "color_shift": color_shift}

mu_best = best["mu"]
ebv_best = best["color_shift"]
d_pc = 10 ** ((mu_best + 5) / 5)
d_kpc = d_pc / 1000

print("\nApproximate NGC 1866 MS fitting result:")
print(f"Best distance modulus mu = {mu_best:.3f} mag")
print(f"Best color shift E(B-V) approx = {ebv_best:.3f} mag")
print(f"Distance = {d_kpc:.2f} kpc")
print("Expected LMC distance ~50 kpc, mu ~18.49 mag")

# -----------------------------
# Overlay shifted MS
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(ngc["B-V"], ngc["Vmag"], s=2, alpha=0.20, label="NGC 1866 stars")
plt.scatter(fit_data["B-V"], fit_data["Vmag"], s=4, alpha=0.35, label="Fit region")

plt.plot(
    hy_binned["B_V"] + ebv_best,
    hy_binned["M_V"] + mu_best,
    linewidth=2.5,
    label=f"Hyades MS shifted: μ={mu_best:.2f}, E(B-V)={ebv_best:.2f}"
)

plt.gca().invert_yaxis()
plt.xlabel("B - V")
plt.ylabel("V")
plt.title("Approximate MS Fitting: NGC 1866")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "ngc1866_ms_fit_overlay.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Compare with Cepheid distance result
# -----------------------------
cepheid_result = ROOT / "cepheid" / "cepheid_distance_result.txt"

print("\nCompare this MS fitting distance with your Cepheid result:")
print("Cepheid result file:", cepheid_result)

out_file = ROOT / "ngc1866_ms" / "ngc1866_ms_fit_result.txt"
with open(out_file, "w") as f:
    f.write("Approximate NGC 1866 MS fitting using Hyades Gaia-calibrated MS\n")
    f.write("WARNING: Approximate because Gaia BP-RP/M_G was transformed roughly to Johnson B-V/V.\n")
    f.write(f"mu = {mu_best:.3f} mag\n")
    f.write(f"E(B-V) approx = {ebv_best:.3f} mag\n")
    f.write(f"Distance = {d_kpc:.2f} kpc\n")
    f.write("Expected LMC distance ~49.9 kpc, mu ~18.49 mag\n")

print("Saved:", out_file)