import numpy as np
import pandas as pd

from config import PHOT_DIR, POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def estimate_snr_and_errors(chosen_aperture=5):
    ensure_directories(OUTPUT_SUBDIRS)

    phot = pd.read_csv(PHOT_DIR / "aperture_photometry.csv")
    pol = pd.read_csv(POL_DIR / "polarimetry_results.csv")

    phot = phot[phot["r_ap"] == chosen_aperture].copy()

    rows = []

    for _, prow in pol.iterrows():
        pid = prow["pair_id"]
        s = phot[phot["pair_id"] == pid].copy()

        # crude mean signal across all beams and HWP states
        signal_vals = []
        sky_vals = []

        for _, row in s.iterrows():
            signal_vals.extend([row["flux_o"], row["flux_e"]])
            sky_vals.extend([row["sky_o"], row["sky_e"]])

        signal_vals = np.array(signal_vals, dtype=float)
        sky_vals = np.array(sky_vals, dtype=float)

        signal_mean = np.nanmean(signal_vals[signal_vals > 0]) if np.any(signal_vals > 0) else np.nan
        sky_mean = np.nanmean(sky_vals)

        if np.isfinite(signal_mean) and signal_mean > 0:
            snr = np.sqrt(signal_mean)
            sigma_p = 1.0 / snr
            sigma_theta_deg = 28.65 * sigma_p / prow["P"] if prow["P"] > 0 else np.nan
        else:
            snr = np.nan
            sigma_p = np.nan
            sigma_theta_deg = np.nan

        rows.append({
            "pair_id": pid,
            "snr_est": snr,
            "sigma_P_est": sigma_p,
            "sigma_theta_deg_est": sigma_theta_deg
        })

    err_df = pd.DataFrame(rows)

    merged = pol.merge(err_df, on="pair_id", how="left")
    out_csv = POL_DIR / "polarimetry_results_with_errors.csv"
    merged.to_csv(out_csv, index=False)

    print(f"Saved: {out_csv}")
    print(merged.head())


if __name__ == "__main__":
    estimate_snr_and_errors()