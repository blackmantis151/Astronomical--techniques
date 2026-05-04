from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

# -----------------------------
# Load Tully-Fisher data
# -----------------------------
tf1 = pd.read_csv(ROOT / "tf" / "J_ApJ_902_145_table1.csv")
tf4 = pd.read_csv(ROOT / "tf" / "J_ApJ_902_145_table4.csv")

df = pd.merge(tf1, tf4, on="PGC")

print("Columns:")
print(df.columns)

# -----------------------------
# Clean and select usable galaxies
# -----------------------------
needed = ["PGC", "Inc", "Wmx", "e_Wmx", "W1mag", "DMbest", "e_DMbest", "Vcmb"]
df = df.dropna(subset=needed).copy()

# Physical cuts
df = df[
    (df["Inc"] > 45) &          # remove face-on galaxies
    (df["Inc"] < 90) &
    (df["Wmx"] > 50) &          # remove very small/invalid linewidths
    (df["e_Wmx"] < 30) &        # keep reasonably measured linewidths
    (df["W1mag"] > 0) &
    (df["DMbest"] > 20) &
    (df["DMbest"] < 40) &
    (df["e_DMbest"] < 0.8)
].copy()

print("\nAfter quality cuts:")
print("Number of galaxies:", len(df))

# -----------------------------
# Tully-Fisher variables
# -----------------------------
# Wmx is corrected HI linewidth. It is roughly 2 Vrot.
# For fitting, catalogs commonly use log linewidth directly.
df["logW"] = np.log10(df["Wmx"])

# Absolute magnitude from catalog distance modulus
df["M_W1"] = df["W1mag"] - df["DMbest"]

# Remove outliers in magnitude space
df = df[
    (df["logW"] > 1.8) &
    (df["logW"] < 3.0) &
    (df["M_W1"] > -26) &
    (df["M_W1"] < -14)
].copy()

print("\nAfter TF range cuts:")
print("Number of galaxies:", len(df))

print("\nSample:")
print(df[["PGC", "Name", "Inc", "Wmx", "logW", "W1mag", "DMbest", "M_W1"]].head())

# -----------------------------
# Fit TF relation
# M = a logW + b
# -----------------------------
coeff = np.polyfit(df["logW"], df["M_W1"], 1)
a, b = coeff

df["M_pred"] = a * df["logW"] + b
df["residual_M"] = df["M_W1"] - df["M_pred"]

# Robust scatter
scatter = 1.4826 * np.median(np.abs(df["residual_M"] - np.median(df["residual_M"])))

print("\nTully-Fisher relation:")
print(f"M_W1 = {a:.3f} log10(Wmx) + {b:.3f}")
print(f"Robust scatter = {scatter:.3f} mag")

# -----------------------------
# Plot 1: TF relation
# -----------------------------
x = np.linspace(df["logW"].min(), df["logW"].max(), 200)
y = a * x + b

plt.figure(figsize=(7, 6))
plt.scatter(df["logW"], df["M_W1"], s=5, alpha=0.45, label="Galaxies")
plt.plot(x, y, linewidth=2.5, label=f"Fit: M={a:.2f} logW {b:+.2f}")

plt.gca().invert_yaxis()
plt.xlabel(r"$\log_{10}(W_{\rm mx})$")
plt.ylabel(r"$M_{W1}$")
plt.title("Tully-Fisher Relation with Inclination Filter")
plt.legend()
plt.grid(True)

plt.savefig(PLOT_DIR / "tully_fisher_relation_inclination_filtered.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 2: residuals
# -----------------------------
plt.figure(figsize=(7, 5))
plt.scatter(df["logW"], df["residual_M"], s=5, alpha=0.45)
plt.axhline(0, linewidth=1)
plt.xlabel(r"$\log_{10}(W_{\rm mx})$")
plt.ylabel(r"$M_{\rm obs} - M_{\rm fit}$")
plt.title("Tully-Fisher Residuals")
plt.grid(True)

plt.savefig(PLOT_DIR / "tully_fisher_residuals.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Recompute TF distances from fitted relation
# -----------------------------
df["mu_tf"] = df["W1mag"] - df["M_pred"]
df["d_tf_mpc"] = 10 ** ((df["mu_tf"] - 25) / 5)

# Catalog distances from DMbest
df["d_catalog_mpc"] = 10 ** ((df["DMbest"] - 25) / 5)

# Remove extreme bad inferred distances
df_dist = df[
    (df["d_tf_mpc"] > 0.1) &
    (df["d_tf_mpc"] < 500) &
    (df["d_catalog_mpc"] > 0.1) &
    (df["d_catalog_mpc"] < 500)
].copy()

print("\nDistance comparison:")
print(df_dist[["d_tf_mpc", "d_catalog_mpc"]].describe())

# -----------------------------
# Plot 3: TF predicted distance vs catalog distance
# -----------------------------
plt.figure(figsize=(6, 6))
plt.scatter(df_dist["d_catalog_mpc"], df_dist["d_tf_mpc"], s=5, alpha=0.45)

lim_min = 0
lim_max = max(df_dist["d_catalog_mpc"].max(), df_dist["d_tf_mpc"].max())
plt.plot([lim_min, lim_max], [lim_min, lim_max], linewidth=2, label="1:1 line")

plt.xlabel("Catalog distance [Mpc]")
plt.ylabel("TF predicted distance [Mpc]")
plt.title("Tully-Fisher Distance Comparison")
plt.legend()
plt.grid(True)

plt.savefig(PLOT_DIR / "tf_distance_comparison_inclination_filtered.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Optional: Hubble-style plot using TF catalog distances
# -----------------------------
# Vcmb is CMB-frame recession velocity in km/s.
hubble = df_dist[
    (df_dist["Vcmb"] > 500) &
    (df_dist["Vcmb"] < 20000) &
    (df_dist["d_catalog_mpc"] > 5)
].copy()

if len(hubble) > 10:
    H0_origin = np.sum(hubble["d_catalog_mpc"] * hubble["Vcmb"]) / np.sum(hubble["d_catalog_mpc"]**2)

    coeff_h = np.polyfit(hubble["d_catalog_mpc"], hubble["Vcmb"], 1)
    H0_fit, intercept = coeff_h

    print("\nHubble-style check using TF catalog distances:")
    print(f"H0 through origin = {H0_origin:.2f} km/s/Mpc")
    print(f"H0 with intercept = {H0_fit:.2f} km/s/Mpc")
    print(f"Intercept = {intercept:.2f} km/s")

    xh = np.linspace(hubble["d_catalog_mpc"].min(), hubble["d_catalog_mpc"].max(), 200)

    plt.figure(figsize=(7, 6))
    plt.scatter(hubble["d_catalog_mpc"], hubble["Vcmb"], s=5, alpha=0.4)
    plt.plot(xh, H0_origin * xh, linewidth=2.5, label=f"Origin fit: H0={H0_origin:.1f}")
    plt.xlabel("TF catalog distance [Mpc]")
    plt.ylabel(r"$V_{\rm CMB}$ [km/s]")
    plt.title("Hubble-style Plot from Tully-Fisher Distances")
    plt.legend()
    plt.grid(True)

    plt.savefig(PLOT_DIR / "tf_hubble_style_plot.png", dpi=300, bbox_inches="tight")
    plt.show()

# -----------------------------
# Save results
# -----------------------------
out = ROOT / "tf" / "tf_results_inclination_filtered.csv"
df_dist.to_csv(out, index=False)

result_file = ROOT / "tf" / "tf_fit_result.txt"
with open(result_file, "w") as f:
    f.write("Tully-Fisher Analysis with Inclination Filter\n")
    f.write("Relation: M_W1 = a log10(Wmx) + b\n")
    f.write(f"Number of galaxies after cuts = {len(df)}\n")
    f.write(f"a = {a:.6f}\n")
    f.write(f"b = {b:.6f}\n")
    f.write(f"Robust scatter = {scatter:.6f} mag\n")
    f.write("\nCuts used:\n")
    f.write("45 < Inc < 90 deg\n")
    f.write("Wmx > 50\n")
    f.write("e_Wmx < 30\n")
    f.write("20 < DMbest < 40\n")
    f.write("e_DMbest < 0.8\n")
    if len(hubble) > 10:
        f.write("\nHubble-style check:\n")
        f.write(f"H0 through origin = {H0_origin:.3f} km/s/Mpc\n")
        f.write(f"H0 with intercept = {H0_fit:.3f} km/s/Mpc\n")
        f.write(f"Intercept = {intercept:.3f} km/s\n")

print("\nSaved:")
print(out)
print(result_file)
print("Plots saved in:", PLOT_DIR)