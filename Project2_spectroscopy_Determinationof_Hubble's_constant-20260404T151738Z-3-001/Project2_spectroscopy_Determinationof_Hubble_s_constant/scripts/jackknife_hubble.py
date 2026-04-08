from pathlib import Path
import json
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
TEXT_DIR = BASE_DIR / "output" / "text"

files = sorted(TEXT_DIR.glob("*_redshift_result.json"))

data = []
for f in files:
    with open(f, "r", encoding="utf-8") as fp:
        d = json.load(fp)
    data.append({
        "filename": d["filename"],
        "distance_mpc": d["distance_mpc"],
        "distance_err": d["distance_err"],
        "velocity_kms": d["velocity_mean_kms"],
        "velocity_err_kms": d["velocity_err_kms"],
    })

def fit_through_origin(distances, velocities, vel_errors):
    distances = np.asarray(distances, dtype=float)
    velocities = np.asarray(velocities, dtype=float)
    vel_errors = np.asarray(vel_errors, dtype=float)

    weights = 1.0 / (vel_errors ** 2)
    H0 = np.sum(weights * distances * velocities) / np.sum(weights * distances**2)
    H0_err = np.sqrt(1.0 / np.sum(weights * distances**2))
    return H0, H0_err

# Full fit
all_d = [x["distance_mpc"] for x in data]
all_v = [x["velocity_kms"] for x in data]
all_ve = [x["velocity_err_kms"] for x in data]

H0_full, H0_full_err = fit_through_origin(all_d, all_v, all_ve)

print("\nFULL SAMPLE")
print("=" * 60)
print(f"H0 = {H0_full:.4f} +/- {H0_full_err:.4f} km/s/Mpc")

# Jackknife
results = []
print("\nJACKKNIFE RESULTS")
print("=" * 60)

for i in range(len(data)):
    subset = data[:i] + data[i+1:]

    d = [x["distance_mpc"] for x in subset]
    v = [x["velocity_kms"] for x in subset]
    ve = [x["velocity_err_kms"] for x in subset]

    H0_i, H0_i_err = fit_through_origin(d, v, ve)
    delta = H0_i - H0_full

    results.append({
        "left_out": data[i]["filename"],
        "H0": H0_i,
        "H0_err": H0_i_err,
        "delta_from_full": delta,
    })

    print(f"Leave out: {data[i]['filename']}")
    print(f"  H0 = {H0_i:.4f} +/- {H0_i_err:.4f}")
    print(f"  Delta = {delta:+.4f}")
    print("-" * 40)

# Optional summary
H0_vals = np.array([r["H0"] for r in results], dtype=float)
print("\nSUMMARY")
print("=" * 60)
print(f"Jackknife mean H0 = {np.mean(H0_vals):.4f}")
print(f"Jackknife std  H0 = {np.std(H0_vals, ddof=1):.4f}")