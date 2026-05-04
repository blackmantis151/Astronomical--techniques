from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

SN_FILE = ROOT / "snia" / "PantheonPlus_Data" / "Pantheon+_Data" / "4_DISTANCES_AND_COVAR" / "Pantheon+SH0ES.dat"

sn = pd.read_csv(SN_FILE, sep=r"\s+", comment="#")

print("Columns:")
print(sn.columns)

print("\nFirst rows:")
print(sn.head())

# Choose distance modulus column
if "MU_SH0ES" in sn.columns:
    mu_col = "MU_SH0ES"
elif "MU" in sn.columns:
    mu_col = "MU"
else:
    raise ValueError("No distance modulus column found. Check printed columns.")

# Choose name column
if "CID" in sn.columns:
    name_col = "CID"
elif "Name" in sn.columns:
    name_col = "Name"
else:
    name_col = None

df = sn.copy()
df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=[mu_col])
df = df[df[mu_col] > 0].copy()

# Distance from distance modulus
df["distance_mpc"] = 10 ** ((df[mu_col] - 25) / 5)

print("\nUsing distance modulus column:", mu_col)
print("\nDistance statistics:")
print(df["distance_mpc"].describe())

# Save cleaned distance table
out_csv = ROOT / "snia" / "snia_distances_from_standard_candle.csv"
df.to_csv(out_csv, index=False)

# Print sample objects
print("\nExample SN Ia distances:")
cols = [mu_col, "distance_mpc"]
if name_col:
    cols = [name_col] + cols

print(df[cols].head(20))

# -----------------------------
# Plot 1: distance modulus distribution
# -----------------------------
plt.figure(figsize=(7,5))
plt.hist(df[mu_col], bins=50)
plt.xlabel("Distance modulus μ")
plt.ylabel("Number of SN Ia")
plt.title("Pantheon+ Type Ia Supernova Distance Moduli")
plt.grid(True)
plt.savefig(PLOT_DIR / "snia_distance_modulus_hist.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 2: distance distribution
# -----------------------------
plt.figure(figsize=(7,5))
plt.hist(df["distance_mpc"], bins=50)
plt.xlabel("Distance [Mpc]")
plt.ylabel("Number of SN Ia")
plt.title("Distances from Type Ia Supernova Standard Candles")
plt.grid(True)
plt.savefig(PLOT_DIR / "snia_distance_hist.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 3: log distance distribution
# -----------------------------
plt.figure(figsize=(7,5))
plt.hist(np.log10(df["distance_mpc"]), bins=50)
plt.xlabel(r"$\log_{10}(d/\mathrm{Mpc})$")
plt.ylabel("Number of SN Ia")
plt.title("Log Distance Distribution of Type Ia Supernovae")
plt.grid(True)
plt.savefig(PLOT_DIR / "snia_log_distance_hist.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Optional: μ vs redshift if redshift exists
# This is not needed for standard candle distance, only context.
# -----------------------------
z_col = None
for c in ["zHD", "zHEL", "zCMB"]:
    if c in df.columns:
        z_col = c
        break

if z_col is not None:
    plt.figure(figsize=(7,6))
    plt.scatter(df[z_col], df[mu_col], s=8, alpha=0.5)
    plt.xscale("log")
    plt.xlabel("Redshift z")
    plt.ylabel("Distance modulus μ")
    plt.title("SN Ia Distance Modulus vs Redshift")
    plt.grid(True, which="both")
    plt.savefig(PLOT_DIR / "snia_mu_vs_redshift_context.png", dpi=300, bbox_inches="tight")
    plt.show()

# Save result summary
result_file = ROOT / "snia" / "snia_standard_candle_distance_result.txt"

with open(result_file, "w") as f:
    f.write("Type Ia Supernova Standard Candle Distance Calculation\n")
    f.write(f"Input file: {SN_FILE}\n")
    f.write(f"Distance modulus column used: {mu_col}\n")
    f.write(f"Number of usable SN Ia: {len(df)}\n")
    f.write("\nDistance statistics in Mpc:\n")
    f.write(str(df["distance_mpc"].describe()))
    f.write("\n\nFormula used:\n")
    f.write("d_Mpc = 10^((mu - 25)/5)\n")

print("\nSaved:")
print(out_csv)
print(result_file)
print("Plots saved in:", PLOT_DIR)


# -----------------------------
# Hubble constant estimation
# -----------------------------
c = 3e5  # km/s

# choose best redshift column
if "zHD" in df.columns:
    z_used = "zHD"
elif "zCMB" in df.columns:
    z_used = "zCMB"
else:
    raise ValueError("No suitable redshift column found")

# Avoid very low z (peculiar velocity dominated)
df_hubble = df[df[z_used] > 0.01].copy()

# velocity
df_hubble["v_kms"] = c * df_hubble[z_used]

# Hubble constant per SN
df_hubble["H0_i"] = df_hubble["v_kms"] / df_hubble["distance_mpc"]

# Clean extreme outliers
df_hubble = df_hubble[(df_hubble["H0_i"] > 30) & (df_hubble["H0_i"] < 100)]

# Final estimate
H0_mean = df_hubble["H0_i"].mean()
H0_std = df_hubble["H0_i"].std()

print("\nEstimated Hubble Constant:")
print(f"H0 = {H0_mean:.2f} ± {H0_std:.2f} km/s/Mpc")

# Save result
with open(result_file, "a") as f:
    f.write("\n\nHubble constant estimation:\n")
    f.write(f"Using redshift column: {z_used}\n")
    f.write(f"Number of SN used: {len(df_hubble)}\n")
    f.write(f"H0 = {H0_mean:.2f} ± {H0_std:.2f} km/s/Mpc\n")

# -----------------------------
# Plot: Hubble diagram (v vs d)
# -----------------------------
plt.figure(figsize=(7,6))
plt.scatter(df_hubble["distance_mpc"], df_hubble["v_kms"], s=8, alpha=0.5)

# Fit line
coeff = np.polyfit(df_hubble["distance_mpc"], df_hubble["v_kms"], 1)
x = np.linspace(df_hubble["distance_mpc"].min(), df_hubble["distance_mpc"].max(), 100)
y = coeff[0]*x + coeff[1]

plt.plot(x, y, color='red', label=f"Fit: H0 ≈ {coeff[0]:.2f}")

plt.xlabel("Distance (Mpc)")
plt.ylabel("Velocity (km/s)")
plt.title("Hubble Diagram from SN Ia")
plt.legend()
plt.grid(True)

plt.savefig(PLOT_DIR / "hubble_diagram_snia.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"\nHubble diagram saved at: {PLOT_DIR / 'hubble_diagram_snia.png'}")