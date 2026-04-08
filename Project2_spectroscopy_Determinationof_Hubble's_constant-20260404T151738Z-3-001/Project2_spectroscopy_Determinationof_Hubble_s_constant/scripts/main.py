from pathlib import Path
import numpy as np

from functions import (
    load_sdss_spectrum,
    get_flux_stats,
    estimate_line_center,
    compute_velocity_from_line,
    predict_observed_wavelength,
    weighted_linear_fit,
)
from plotting import (
    plot_full_spectrum,
    plot_line_region,
    plot_hubble_diagram,
)

# --------------------------------------------------
# Paths
# --------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_hubble"
PLOT_DIR = BASE_DIR / "output" / "plots"


# --------------------------------------------------
# Galaxy information from the manual
# Fill these carefully from your sheet
# filename : distance, distance_err, reference_line_name, rest_lambda
# --------------------------------------------------
galaxies = [
    {
        "filename": "spec-0501-52235-0386.fits",
        "distance_mpc": None,
        "distance_err": None,
        "ref_line_name": "Halpha",
        "ref_rest_lambda": 6562.8,
        "line_mode": "emission",
    },
    # Add the remaining 7 galaxies here
]


# --------------------------------------------------
# Candidate rest wavelengths to confirm the redshift
# Add or remove as needed
# --------------------------------------------------
candidate_lines = {
    "Halpha": 6562.8,
    "Hbeta": 4861.3,
    "OIII_4959": 4958.9,
    "OIII_5007": 5006.8,
    "NII_6548": 6548.1,
    "NII_6583": 6583.6,
    "SII_6716": 6716.4,
    "SII_6731": 6730.8,
    "MgII": 2798.0,
    "CIV": 1549.0,
}


def analyze_one_galaxy(gal):
    """
    Analyze one galaxy using its reference line first.
    """
    filepath = DATA_DIR / gal["filename"]
    wavelength, flux, header = load_sdss_spectrum(filepath)

    print("\n" + "=" * 60)
    print(f"File: {gal['filename']}")
    print(f"Distance (Mpc): {gal['distance_mpc']}")
    print(f"Reference line: {gal['ref_line_name']}  Rest lambda = {gal['ref_rest_lambda']} A")

    stats = get_flux_stats(flux)
    print("Flux stats:", stats)

    # Step 1: plot full spectrum
    plot_full_spectrum(
        wavelength,
        flux,
        title=gal["filename"],
        savepath=PLOT_DIR / f"{gal['filename']}_full.png",
        show=True,
    )

    # Step 2: rough detection of reference line
    obs_lambda, obs_flux, local_data = estimate_line_center(
        wavelength,
        flux,
        expected_center=gal["ref_rest_lambda"],   # first rough try
        window=200,
        mode=gal["line_mode"],
    )

    # If this spectrum is redshifted a lot, this may fail.
    # Then you should inspect the plot and manually adjust the expected center.
    if obs_lambda is None:
        print("Could not detect the reference line.")
        return None

    z_ref, v_ref = compute_velocity_from_line(obs_lambda, gal["ref_rest_lambda"])

    print(f"Detected reference line at: {obs_lambda:.2f} A")
    print(f"Estimated z from reference line: {z_ref:.5f}")
    print(f"Estimated velocity from reference line: {v_ref:.2f} km/s")

    w_cut, f_cut = local_data
    plot_line_region(
        w_cut,
        f_cut,
        expected_center=gal["ref_rest_lambda"],
        detected_center=obs_lambda,
        title=f"{gal['filename']} : {gal['ref_line_name']}",
        savepath=PLOT_DIR / f"{gal['filename']}_{gal['ref_line_name']}.png",
        show=True,
    )

    # Step 3: confirm using other lines
    measured_velocities = [v_ref]

    for line_name, rest_lambda in candidate_lines.items():
        if line_name == gal["ref_line_name"]:
            continue

        predicted = predict_observed_wavelength(rest_lambda, z_ref)

        obs2, flux2, local2 = estimate_line_center(
            wavelength,
            flux,
            expected_center=predicted,
            window=30,
            mode=gal["line_mode"],
        )

        if obs2 is None:
            continue

        z2, v2 = compute_velocity_from_line(obs2, rest_lambda)

        # keep only roughly consistent lines
        if abs(v2 - v_ref) < 3000:
            measured_velocities.append(v2)

            w2, f2 = local2
            plot_line_region(
                w2,
                f2,
                expected_center=predicted,
                detected_center=obs2,
                title=f"{gal['filename']} : {line_name}",
                savepath=PLOT_DIR / f"{gal['filename']}_{line_name}.png",
                show=False,
            )

            print(f"Confirmed line: {line_name:10s}  obs = {obs2:.2f} A  v = {v2:.2f} km/s")

    measured_velocities = np.array(measured_velocities)
    v_mean = np.mean(measured_velocities)

    if len(measured_velocities) > 1:
        v_err = np.std(measured_velocities, ddof=1)
    else:
        v_err = 500.0   # fallback rough uncertainty

    print(f"Mean velocity = {v_mean:.2f} km/s")
    print(f"Velocity error = {v_err:.2f} km/s")
    print(f"Number of accepted lines = {len(measured_velocities)}")

    return {
        "filename": gal["filename"],
        "distance_mpc": gal["distance_mpc"],
        "distance_err": gal["distance_err"],
        "velocity_kms": v_mean,
        "velocity_err": v_err,
    }


def main():
    results = []

    for gal in galaxies:
        result = analyze_one_galaxy(gal)
        if result is not None:
            results.append(result)

    # keep only galaxies with known distances
    results = [r for r in results if r["distance_mpc"] is not None]

    if len(results) < 2:
        print("\nNot enough galaxies with distances filled in for Hubble fit.")
        return

    distance = np.array([r["distance_mpc"] for r in results], dtype=float)
    velocity = np.array([r["velocity_kms"] for r in results], dtype=float)
    velocity_err = np.array([r["velocity_err"] for r in results], dtype=float)

    m, b, m_err, b_err = weighted_linear_fit(distance, velocity, velocity_err)

    print("\n" + "=" * 60)
    print("Weighted Hubble fit")
    print(f"H0 = {m:.2f} +/- {m_err:.2f} km/s/Mpc")
    print(f"Intercept = {b:.2f} +/- {b_err:.2f} km/s")

    # Age estimate in Gyr
    age_gyr = 978.0 / m
    age_err_gyr = 978.0 * m_err / (m ** 2)

    print(f"Age of universe = {age_gyr:.2f} +/- {age_err_gyr:.2f} Gyr")

    plot_hubble_diagram(
        distance,
        velocity,
        velocity_err=velocity_err,
        slope=m,
        intercept=b,
        title="Hubble Diagram",
        savepath=PLOT_DIR / "hubble_diagram.png",
        show=True,
    )


if __name__ == "__main__":
    main()