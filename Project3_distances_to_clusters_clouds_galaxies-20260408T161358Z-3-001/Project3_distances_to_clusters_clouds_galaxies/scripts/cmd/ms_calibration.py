from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

# Load selected Hyades members
df = pd.read_csv(ROOT / "gaia" / "hyades_members_selected.csv")

# Clean
df = df.dropna(subset=["parallax", "phot_g_mean_mag", "bp_rp"])
df = df[df["parallax"] > 0].copy()

# Distance from Gaia parallax
df["distance_pc"] = 1000.0 / df["parallax"]

# Absolute Gaia G magnitude
df["M_G"] = df["phot_g_mean_mag"] - 5*np.log10(df["distance_pc"]) + 5

print("Hyades absolute CMD data:")
print(df[["bp_rp", "phot_g_mean_mag", "parallax", "distance_pc", "M_G"]].head())

print("\nDistance statistics from selected Hyades members:")
print(df["distance_pc"].describe())

print("\nMedian Hyades distance:")
print(np.median(df["distance_pc"]), "pc")

# Save table with absolute magnitudes
out_csv = ROOT / "gaia" / "hyades_members_absolute_cmd.csv"
df.to_csv(out_csv, index=False)
print(f"\nSaved absolute CMD table: {out_csv}")

# Plot apparent CMD
plt.figure(figsize=(7,6))
plt.scatter(df["bp_rp"], df["phot_g_mean_mag"], s=12, alpha=0.85)
plt.gca().invert_yaxis()
plt.xlabel("BP - RP")
plt.ylabel("G")
plt.title("Hyades Apparent CMD")
plt.grid(True)
plt.savefig(PLOT_DIR / "hyades_apparent_cmd_selected.png", dpi=300, bbox_inches="tight")
plt.show()

# Plot absolute CMD
plt.figure(figsize=(7,6))
plt.scatter(df["bp_rp"], df["M_G"], s=12, alpha=0.85)
plt.gca().invert_yaxis()
plt.xlabel("BP - RP")
plt.ylabel(r"$M_G$")
plt.title("Hyades Absolute CMD from Gaia Parallax")
plt.grid(True)
plt.savefig(PLOT_DIR / "hyades_absolute_cmd_parallax.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------------------
# Overlay Apparent vs Absolute CMD
# -----------------------------------------

plt.figure(figsize=(7,6))

# Apparent CMD
plt.scatter(df["bp_rp"], df["phot_g_mean_mag"],
            s=10, alpha=0.5, label="Apparent (G)")

# Absolute CMD
plt.scatter(df["bp_rp"], df["M_G"],
            s=10, alpha=0.8, label="Absolute ($M_G$)")

plt.gca().invert_yaxis()
plt.xlabel("BP - RP")
plt.ylabel("Magnitude")
plt.title("Hyades CMD: Apparent vs Absolute")

plt.legend()
plt.grid(True)

plt.savefig(PLOT_DIR / "hyades_cmd_overlay.png", dpi=300, bbox_inches="tight")
plt.show()