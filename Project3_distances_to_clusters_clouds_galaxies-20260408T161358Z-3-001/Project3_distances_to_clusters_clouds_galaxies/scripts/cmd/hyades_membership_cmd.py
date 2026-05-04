from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

# Load Hyades Gaia data
df = pd.read_csv(ROOT / "gaia" / "gaia_hyades.csv")

# Basic clean
df = df.dropna(subset=["ra", "dec", "parallax", "pmra", "pmdec", "phot_g_mean_mag", "bp_rp", "ruwe"])
df = df[(df["parallax"] > 0) & (df["ruwe"] < 1.4)].copy()

# -----------------------------
# 1. Vector Point Diagram
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(df["pmra"], df["pmdec"], s=5, alpha=0.35)
plt.xlabel(r"$\mu_{\alpha *}$ = pmRA [mas/yr]")
plt.ylabel(r"$\mu_\delta$ = pmDec [mas/yr]")
plt.title("Hyades: Vector Point Diagram")
plt.grid(True)
plt.savefig(PLOT_DIR / "hyades_vpd_all.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 2. Select likely Hyades members
# Approx Hyades mean motion:
# pmra ~ +100 mas/yr, pmdec ~ -25 mas/yr
# parallax ~ 21 mas
# These are broad cuts; adjust after seeing the plot.
# -----------------------------
members = df[
    (df["pmra"] > 50) & (df["pmra"] < 160) &
    (df["pmdec"] > -80) & (df["pmdec"] < 30) &
    (df["parallax"] > 15) & (df["parallax"] < 30)
].copy()

field = df.loc[~df.index.isin(members.index)].copy()

print("Total stars:", len(df))
print("Likely Hyades members:", len(members))
print("Field stars:", len(field))

print("\nMember parallax stats:")
print(members["parallax"].describe())

print("\nMember distance estimate:")
members["distance_pc"] = 1000 / members["parallax"]
print(members["distance_pc"].describe())

print("\nMedian distance from parallax:")
print(np.median(members["distance_pc"]), "pc")

# Save selected members
members.to_csv(ROOT / "gaia" / "hyades_members_selected.csv", index=False)

# -----------------------------
# 3. VPD with selected members
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(field["pmra"], field["pmdec"], s=4, alpha=0.2, label="Field")
plt.scatter(members["pmra"], members["pmdec"], s=10, alpha=0.9, label="Selected Hyades")
plt.xlabel(r"$\mu_{\alpha *}$ = pmRA [mas/yr]")
plt.ylabel(r"$\mu_\delta$ = pmDec [mas/yr]")
plt.title("Hyades: VPD Member Selection")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "hyades_vpd_members.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 4. Sky plot
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(field["ra"], field["dec"], s=4, alpha=0.15, label="Field")
plt.scatter(members["ra"], members["dec"], s=10, alpha=0.9, label="Selected Hyades")
plt.xlabel("RA [deg]")
plt.ylabel("Dec [deg]")
plt.title("Hyades: Sky Distribution")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "hyades_sky_members.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 5. CMD of all stars vs selected members
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(field["bp_rp"], field["phot_g_mean_mag"], s=4, alpha=0.15, label="Field")
plt.scatter(members["bp_rp"], members["phot_g_mean_mag"], s=12, alpha=0.9, label="Selected Hyades")
plt.gca().invert_yaxis()
plt.xlabel("BP - RP")
plt.ylabel("G")
plt.title("Hyades CMD after Proper Motion + Parallax Selection")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "hyades_cmd_members.png", dpi=300, bbox_inches="tight")
plt.show()