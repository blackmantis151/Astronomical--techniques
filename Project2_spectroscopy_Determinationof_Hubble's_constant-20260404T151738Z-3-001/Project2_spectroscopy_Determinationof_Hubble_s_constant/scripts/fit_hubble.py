from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt

from functions import weighted_linear_fit

BASE_DIR = Path(__file__).resolve().parent.parent
TEXT_DIR = BASE_DIR / "output" / "text"
PLOT_DIR = BASE_DIR / "output" / "plots"

PLOT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# Exclude problematic objects
# --------------------------------------------------
EXCLUDE_FILES = {}

# --------------------------------------------------
# Step 1: load all results
# --------------------------------------------------
files = list(TEXT_DIR.glob("*_redshift_result.json"))

distances = []
dist_errors = []
velocities = []
vel_errors = []
used_files = []

for f in files:
    with open(f, "r", encoding="utf-8") as fp:
        data = json.load(fp)

    filename = data["filename"]

    if filename in EXCLUDE_FILES:
        continue

    d = data["distance_mpc"]
    d_err = data["distance_err"]
    v = data["velocity_mean_kms"]
    v_err = data["velocity_err_kms"]

    distances.append(d)
    dist_errors.append(d_err)
    velocities.append(v)
    vel_errors.append(v_err)
    used_files.append(filename)

distances = np.array(distances, dtype=float)
dist_errors = np.array(dist_errors, dtype=float)
velocities = np.array(velocities, dtype=float)
vel_errors = np.array(vel_errors, dtype=float)

print("\nIncluded files:")
for name in used_files:
    print(" ", name)

# --------------------------------------------------
# Step 2: weighted fit (with intercept)
# --------------------------------------------------
m, b, m_err, b_err = weighted_linear_fit(distances, velocities, vel_errors)

print("\nFit results (with intercept):")
print("H0 =", m, "+/-", m_err, "km/s/Mpc")
print("Intercept =", b, "+/-", b_err)

# --------------------------------------------------
# Step 3: force intercept = 0 fit
# --------------------------------------------------
weights = 1.0 / vel_errors**2
H0 = np.sum(weights * distances * velocities) / np.sum(weights * distances**2)
H0_err = np.sqrt(1.0 / np.sum(weights * distances**2))

print("\nFit results (forced through origin):")
print("H0 =", H0, "+/-", H0_err, "km/s/Mpc")

# --------------------------------------------------
# Step 4: plot
# --------------------------------------------------
plt.figure(figsize=(8, 6))

plt.errorbar(
    distances,
    velocities,
    yerr=vel_errors,
    xerr=dist_errors,
    fmt='o',
    capsize=3,
    label="Data"
)

x_fit = np.linspace(0, max(distances) * 1.1, 100)
y_fit = H0 * x_fit
plt.plot(x_fit, y_fit, label=f"H0 = {H0:.2f}")

for x, y, name in zip(distances, velocities, used_files):
    plt.annotate(name.replace(".fits", ""), (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)

plt.xlabel("Distance (Mpc)")
plt.ylabel("Velocity (km/s)")
plt.title("Cleaned Hubble Diagram")
plt.grid(alpha=0.3)
plt.legend()

plt.savefig(PLOT_DIR / "hubble_fit_cleaned.png", dpi=200, bbox_inches="tight")
plt.show()

# --------------------------------------------------
# Step 5: Age of Universe
# --------------------------------------------------
H0_SI = H0 * 1000 / (3.086e22)
age_sec = 1 / H0_SI
age_yr = age_sec / (60 * 60 * 24 * 365.25)

print("\nEstimated age of Universe:")
print("Age =", age_yr / 1e9, "Gyr")