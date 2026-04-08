from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent
TEXT_DIR = BASE_DIR / "output" / "text"
PLOT_DIR = BASE_DIR / "output" / "plots"

PLOT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# Files to exclude from the low-z linear Hubble fit
# --------------------------------------------------
EXCLUDE_FILES = {
    "spec-5193-56066-0534.fits",
    "spec-1617-53112-0467.fits",
    "spec-1818-54539-0147.fits",
    "spec-0501-52235-0386.fits",
}

# --------------------------------------------------
# Load results
# --------------------------------------------------
files = sorted(TEXT_DIR.glob("*_redshift_result.json"))

data = []
for f in files:
    with open(f, "r", encoding="utf-8") as fp:
        d = json.load(fp)

    filename = d["filename"]

    if filename in EXCLUDE_FILES:
        continue

    data.append({
        "filename": filename,
        "distance_mpc": d["distance_mpc"],
        "distance_err": d["distance_err"],
        "velocity_kms": d["velocity_mean_kms"],
        "velocity_err_kms": d["velocity_err_kms"],
    })

if len(data) < 2:
    raise ValueError("Not enough low-z galaxies left for fitting.")

distances = np.array([x["distance_mpc"] for x in data], dtype=float)
dist_errors = np.array([x["distance_err"] for x in data], dtype=float)
velocities = np.array([x["velocity_kms"] for x in data], dtype=float)
vel_errors = np.array([x["velocity_err_kms"] for x in data], dtype=float)

# --------------------------------------------------
# Weighted fit through origin: v = H0 d
# --------------------------------------------------
weights = 1.0 / (vel_errors ** 2)

H0 = np.sum(weights * distances * velocities) / np.sum(weights * distances**2)
H0_formal_err = np.sqrt(1.0 / np.sum(weights * distances**2))

# --------------------------------------------------
# Jackknife uncertainty on low-z sample
# --------------------------------------------------
jackknife_vals = []

for i in range(len(data)):
    d_sub = np.delete(distances, i)
    v_sub = np.delete(velocities, i)
    ve_sub = np.delete(vel_errors, i)

    w_sub = 1.0 / (ve_sub ** 2)
    H0_i = np.sum(w_sub * d_sub * v_sub) / np.sum(w_sub * d_sub**2)
    jackknife_vals.append(H0_i)

jackknife_vals = np.array(jackknife_vals, dtype=float)
H0_jackknife_std = np.std(jackknife_vals, ddof=1) if len(jackknife_vals) > 1 else np.nan

print("\nLOW-z Hubble fit")
print("=" * 60)
print("Included galaxies:")
for x in data:
    print(" ", x["filename"])

print(f"\nH0 (through origin) = {H0:.4f} km/s/Mpc")
print(f"Formal fit error    = {H0_formal_err:.4f} km/s/Mpc")
print(f"Jackknife std       = {H0_jackknife_std:.4f} km/s/Mpc")

# --------------------------------------------------
# Age of Universe
# --------------------------------------------------
# H0 [km/s/Mpc] -> s^-1
H0_SI = H0 * 1000.0 / (3.085677581e22)

age_sec = 1.0 / H0_SI
age_yr = age_sec / (365.25 * 24 * 3600)
age_gyr = age_yr / 1e9

print(f"Age of Universe     = {age_gyr:.4f} Gyr")

# --------------------------------------------------
# Plot
# --------------------------------------------------
plt.figure(figsize=(8, 6))

plt.errorbar(
    distances,
    velocities,
    yerr=vel_errors,
    xerr=dist_errors,
    fmt="o",
    capsize=4,
    label="Low-z galaxies"
)

xfit = np.linspace(0, max(distances) * 1.1, 200)
yfit = H0 * xfit

plt.plot(xfit, yfit, label=f"H0 = {H0:.2f} km/s/Mpc")

for x in data:
    plt.annotate(
        x["filename"].replace(".fits", ""),
        (x["distance_mpc"], x["velocity_kms"]),
        textcoords="offset points",
        xytext=(5, 5),
        fontsize=8
    )

plt.xlabel("Distance (Mpc)")
plt.ylabel("Velocity (km/s)")
plt.title("Low-z Hubble Diagram")
plt.grid(alpha=0.3)
plt.legend()

savepath = PLOT_DIR / "hubble_fit_lowz.png"
plt.savefig(savepath, dpi=200, bbox_inches="tight")
plt.show()

print(f"\nSaved plot: {savepath}")