import numpy as np
import pandas as pd

from config import POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def select_good_polarimetry(min_p_over_sigma=0.8, min_snr=5.0):
    ensure_directories(OUTPUT_SUBDIRS)

    df = pd.read_csv(POL_DIR / "polarimetry_results_with_errors.csv")

    df["p_over_sigma"] = df["P"] / df["sigma_P_est"]

    good = df[
        np.isfinite(df["P"]) &
        np.isfinite(df["theta_deg"]) &
        np.isfinite(df["sigma_P_est"]) &
        np.isfinite(df["snr_est"]) &
        (df["snr_est"] >= min_snr) &
        (df["p_over_sigma"] >= min_p_over_sigma)
    ].copy()

    rejected = df.drop(good.index).copy()

    good_csv = POL_DIR / "polarimetry_results_good.csv"
    rej_csv = POL_DIR / "polarimetry_results_low_quality.csv"

    good.to_csv(good_csv, index=False)
    rejected.to_csv(rej_csv, index=False)

    print(f"Saved good sample: {good_csv}")
    print(f"Saved low-quality sample: {rej_csv}")
    print(f"Good stars: {len(good)}")
    print(f"Low-quality stars: {len(rejected)}")


if __name__ == "__main__":
    select_good_polarimetry()