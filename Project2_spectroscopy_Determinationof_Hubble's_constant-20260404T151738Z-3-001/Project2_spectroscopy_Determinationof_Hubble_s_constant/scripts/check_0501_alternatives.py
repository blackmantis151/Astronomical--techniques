from pathlib import Path
import numpy as np

from functions import (
    load_sdss_spectrum,
    estimate_line_center,
    refine_line_center_centroid,
    compute_velocity_from_line,
    predict_observed_wavelength,
)
from config import REST_LINES

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_hubble"

filename = "spec-0501-52235-0386.fits"
filepath = DATA_DIR / filename

wavelength, flux, header = load_sdss_spectrum(filepath)

# Try alternative assignments for the suspicious features
# format: (label, rest_lambda, guessed_obs_lambda)
tests = [
    ("Halpha_at_6745", REST_LINES["Halpha"], 6745.0),
    ("Halpha_at_6810", REST_LINES["Halpha"], 6810.0),
    ("NII_6544_at_6745", REST_LINES["NII_6544"], 6745.0),
    ("NII_6583_at_6745", REST_LINES["NII_6583"], 6745.0),
]

# Lines to check once a redshift is obtained
check_lines = ["Hbeta", "OIII_5007", "SII_6733"]

for test_name, rest_lambda, guess_obs in tests:
    print("\n" + "=" * 70)
    print("TEST:", test_name)
    print("Rest wavelength:", rest_lambda)
    print("Manual guess:", guess_obs)

    obs_rough, _, _ = estimate_line_center(
        wavelength, flux,
        expected_center=guess_obs,
        window=20,
        mode="emission"
    )

    if obs_rough is None:
        print("No feature found near this guess.")
        continue

    obs_refined, strength, _ = refine_line_center_centroid(
        wavelength, flux,
        expected_center=obs_rough,
        window=15,
        mode="emission"
    )

    if obs_refined is None:
        print("Could not refine feature.")
        continue

    z_ref, v_ref = compute_velocity_from_line(obs_refined, rest_lambda)

    print(f"Refined observed wavelength: {obs_refined:.3f}")
    print(f"z = {z_ref:.6f}")
    print(f"v = {v_ref:.2f} km/s")

    accepted = []
    accepted.append((test_name, obs_refined, z_ref, v_ref))

    for line_name in check_lines:
        rest2 = REST_LINES[line_name]
        pred2 = predict_observed_wavelength(rest2, z_ref)

        obs2_rough, _, _ = estimate_line_center(
            wavelength, flux,
            expected_center=pred2,
            window=20,
            mode="emission"
        )

        if obs2_rough is None:
            print(f"{line_name}: no feature near predicted {pred2:.2f}")
            continue

        obs2_refined, strength2, _ = refine_line_center_centroid(
            wavelength, flux,
            expected_center=obs2_rough,
            window=15,
            mode="emission"
        )

        if obs2_refined is None:
            print(f"{line_name}: refinement failed")
            continue

        z2, v2 = compute_velocity_from_line(obs2_refined, rest2)
        dz = abs(z2 - z_ref)

        print(
            f"{line_name}: predicted {pred2:.2f}, detected {obs2_refined:.2f}, "
            f"z={z2:.6f}, dz={dz:.6f}"
        )

        if dz < 0.003:
            accepted.append((line_name, obs2_refined, z2, v2))

    if len(accepted) > 1:
        velocities = np.array([x[3] for x in accepted], dtype=float)
        print(f"Accepted lines: {[x[0] for x in accepted]}")
        print(f"Mean velocity = {np.mean(velocities):.2f} km/s")
        if len(velocities) > 1:
            print(f"Velocity scatter = {np.std(velocities, ddof=1):.2f} km/s")
    else:
        print("No convincing multi-line confirmation.")