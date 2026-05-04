from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

FJ_FILE = ROOT / "fj" / "J_MNRAS_443_1231_FPsample.csv"

df = pd.read_csv(FJ_FILE)

print("Columns:")
print(df.columns)
print(df.head())

# -----------------------------
# Use K-band total magnitude and velocity dispersion
# -----------------------------
needed = ["Ktot", "logVd", "cz"]
df = df.dropna(subset=needed).copy()

# Adopt H0 from previous ladder/Hubble result.
# You can change this to your SN value 64.54 or TF value 70.51.
H0_CAL = 70.51  # km/s/Mpc

# Clean sample
df = df[
    (df["cz"] > 1000) &
    (df["cz"] < 20000) &
    (df["logVd"] > 1.7) &
    (df["logVd"] < 2.6) &
    (df["Ktot"] > 5) &
    (df["Ktot"] < 16)
].copy()

# Distance from Hubble flow, calibrated by previous ladder H0
df["d_mpc"] = df["cz"] / H0_CAL

# Absolute K magnitude
df["M_K"] = df["Ktot"] - 5*np.log10(df["d_mpc"]) - 25

# Clean absolute magnitude range
df = df[(df["M_K"] > -28) & (df["M_K"] < -18)].copy()

print("\nUsable galaxies:", len(df))
print(df[["Ktot", "logVd", "cz", "d_mpc", "M_K"]].head())

# -----------------------------
# Fit Faber-Jackson relation
# M_K = a log(sigma) + b
# -----------------------------
a, b = np.polyfit(df["logVd"], df["M_K"], 1)

df["M_pred"] = a*df["logVd"] + b
df["residual"] = df["M_K"] - df["M_pred"]

scatter = 1.4826 * np.median(np.abs(df["residual"] - np.median(df["residual"])))

print("\nFaber-Jackson relation:")
print(f"M_K = {a:.3f} log10(sigma) + {b:.3f}")
print(f"Robust scatter = {scatter:.3f} mag")

# -----------------------------
# Plot 1: Faber-Jackson relation
# -----------------------------
x = np.linspace(df["logVd"].min(), df["logVd"].max(), 200)
y = a*x + b

plt.figure(figsize=(7, 6))
plt.scatter(df["logVd"], df["M_K"], s=5, alpha=0.45, label="6dFGS early type galaxies")
plt.plot(x, y, linewidth=2.5, label=f"Fit: M={a:.2f} logσ {b:+.2f}")

plt.gca().invert_yaxis()
plt.xlabel(r"$\log_{10}(\sigma)$")
plt.ylabel(r"$M_K$")
plt.title("Faber-Jackson Relation")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "faber_jackson_relation.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 2: residuals
# -----------------------------
plt.figure(figsize=(7, 5))
plt.scatter(df["logVd"], df["residual"], s=5, alpha=0.45)
plt.axhline(0, linewidth=1)
plt.xlabel(r"$\log_{10}(\sigma)$")
plt.ylabel(r"$M_{\rm obs}-M_{\rm fit}$")
plt.title("Faber-Jackson Residuals")
plt.grid(True)
plt.savefig(PLOT_DIR / "faber_jackson_residuals.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Estimate distances using fitted FJ relation
# -----------------------------
df["mu_fj"] = df["Ktot"] - df["M_pred"]
df["d_fj_mpc"] = 10 ** ((df["mu_fj"] - 25)/5)

# Compare with Hubble-flow distances
plt.figure(figsize=(6, 6))
plt.scatter(df["d_mpc"], df["d_fj_mpc"], s=5, alpha=0.4)

lim = max(df["d_mpc"].max(), df["d_fj_mpc"].max())
plt.plot([0, lim], [0, lim], linewidth=2, label="1:1 line")

plt.xlabel("Hubble-flow distance [Mpc]")
plt.ylabel("Faber-Jackson distance [Mpc]")
plt.title("Faber-Jackson Distance Comparison")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "faber_jackson_distance_comparison.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Save results
# -----------------------------
out_csv = ROOT / "fj" / "fj_results.csv"
df.to_csv(out_csv, index=False)

result_file = ROOT / "fj" / "fj_fit_result.txt"
with open(result_file, "w") as f:
    f.write("Faber-Jackson Analysis\n")
    f.write(f"H0 calibration used = {H0_CAL:.3f} km/s/Mpc\n")
    f.write("Relation: M_K = a log10(sigma) + b\n")
    f.write(f"Number of galaxies = {len(df)}\n")
    f.write(f"a = {a:.6f}\n")
    f.write(f"b = {b:.6f}\n")
    f.write(f"Robust scatter = {scatter:.6f} mag\n")
    f.write("\nFormulae:\n")
    f.write("d_Mpc = cz / H0\n")
    f.write("M_K = Ktot - 5 log10(d_Mpc) - 25\n")
    f.write("mu_FJ = Ktot - M_pred\n")
    f.write("d_FJ = 10^((mu_FJ - 25)/5)\n")

print("\nSaved:")
print(out_csv)
print(result_file)
print("Plots saved in:", PLOT_DIR)