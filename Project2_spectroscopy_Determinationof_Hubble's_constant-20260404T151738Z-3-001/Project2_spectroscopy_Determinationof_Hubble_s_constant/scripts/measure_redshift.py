from pathlib import Path
import json
import numpy as np

from functions import (
    load_sdss_spectrum,
    estimate_line_center,
    refine_line_center_centroid,
    compute_velocity_from_line,
    predict_observed_wavelength,
)
from config import GALAXIES, REST_LINES

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_hubble"
OUT_DIR = BASE_DIR / "output"
TEXT_DIR = OUT_DIR / "text"

TEXT_DIR.mkdir(parents=True, exist_ok=True)

# choose one galaxy for now
INDEX = 0

gal = GALAXIES[INDEX]

filename = gal["filename"]
ref_name = gal["reference_line"]
rest_lambda = gal["rest_lambda"]
guess_obs_lambda = gal["guess_obs_lambda"]
line_mode = gal["line_mode"]

if guess_obs_lambda is None:
    raise ValueError("Please fill guess_obs_lambda in config.py first.")

filepath = DATA_DIR / filename
wavelength, flux, header = load_sdss_spectrum(filepath)

# rough local detection
obs_lambda_rough, obs_flux_rough, local_data = estimate_line_center(
    wavelength,
    flux,
    expected_center=guess_obs_lambda,
    window=25,
    mode=line_mode,
)

# refined center
obs_lambda_refined, feature_strength, local_refined = refine_line_center_centroid(
    wavelength,
    flux,
    expected_center=obs_lambda_rough,
    window=20,
    mode=line_mode,
)

if obs_lambda_refined is None:
    raise ValueError("Could not refine the reference line center.")

z_ref, v_ref = compute_velocity_from_line(obs_lambda_refined, rest_lambda)

print("Reference line result")
print("File:", filename)
print("Reference line:", ref_name)
print("Rest wavelength:", rest_lambda)
print("Manual guess:", guess_obs_lambda)
print("Rough center:", obs_lambda_rough)
print("Refined center:", obs_lambda_refined)
print("Redshift z:", z_ref)
print("Velocity (km/s):", v_ref)

# choose a few likely confirmation lines
if ref_name == "Halpha":
    candidates = ["Hbeta", "OIII_5007", "SII_6733"]
elif ref_name == "CIV":
    candidates = ["CIII", "MgII"]
elif ref_name == "MgII":
    candidates = ["OII", "Hgamma", "Hbeta"]
else:
    candidates = []

results = []
results.append({
    "line_name": ref_name,
    "rest_lambda": rest_lambda,
    "obs_lambda": obs_lambda_refined,
    "z": z_ref,
    "velocity_kms": v_ref,
})

for line_name in candidates:
    rest2 = REST_LINES[line_name]
    pred2 = predict_observed_wavelength(rest2, z_ref)

    obs2_rough, _, _ = estimate_line_center(
        wavelength,
        flux,
        expected_center=pred2,
        window=20,
        mode=line_mode,
    )

    if obs2_rough is None:
        continue

    obs2_refined, strength2, _ = refine_line_center_centroid(
        wavelength,
        flux,
        expected_center=obs2_rough,
        window=15,
        mode=line_mode,
    )

    if obs2_refined is None:
        continue

    z2, v2 = compute_velocity_from_line(obs2_refined, rest2)

    # keep only lines consistent with reference z
    if abs(z2 - z_ref) < 0.003:
        results.append({
            "line_name": line_name,
            "rest_lambda": rest2,
            "obs_lambda": obs2_refined,
            "z": z2,
            "velocity_kms": v2,
        })

velocities = np.array([r["velocity_kms"] for r in results], dtype=float)
v_mean = float(np.mean(velocities))
v_err = float(np.std(velocities, ddof=1)) if len(velocities) > 1 else 500.0

summary = {
    "filename": filename,
    "reference_line": ref_name,
    "distance_mpc": gal["distance_mpc"],
    "distance_err": gal["distance_err"],
    "z_reference": z_ref,
    "velocity_reference_kms": v_ref,
    "accepted_lines": results,
    "velocity_mean_kms": v_mean,
    "velocity_err_kms": v_err,
}

outpath = TEXT_DIR / f"{filename}_redshift_result.json"
with open(outpath, "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

print("\nAccepted lines:")
for r in results:
    print(r)

print("\nMean velocity:", v_mean)
print("Velocity error:", v_err)
print("Saved:", outpath)