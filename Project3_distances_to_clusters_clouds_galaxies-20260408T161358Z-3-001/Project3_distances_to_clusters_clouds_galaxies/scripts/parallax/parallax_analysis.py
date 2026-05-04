from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Project root
ROOT = Path(__file__).resolve().parents[2]

# Load data
df = pd.read_csv(ROOT / "gaia" / "gaia_nearby_ms.csv")

# Clean
df = df[df["parallax"] > 0].copy()

# Distance (pc)
df["distance_pc"] = 1000 / df["parallax"]

# Absolute magnitude
df["M_G"] = df["phot_g_mean_mag"] - 5*np.log10(df["distance_pc"]) + 5

print("\nSample data:")
print(df[["parallax", "distance_pc", "phot_g_mean_mag", "M_G"]].head())

print("\nDistance stats:")
print(df["distance_pc"].describe())

# Plot
plt.figure(figsize=(7,5))
plt.hist(df["distance_pc"], bins=50)
plt.xlabel("Distance (pc)")
plt.ylabel("Number of stars")
plt.title("Distance Distribution from Gaia Parallax")
plt.grid(True)

# Save plot
plot_path = ROOT / "scripts" / "plots" / "parallax_distance_hist.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")

print(f"\nPlot saved at: {plot_path}")

plt.show()