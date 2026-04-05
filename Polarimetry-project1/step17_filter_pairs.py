import pandas as pd
import numpy as np

from config import POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def filter_pairs():
    ensure_directories(OUTPUT_SUBDIRS)

    pair_path = POL_DIR / "beam_pairs.csv"
    df = pd.read_csv(pair_path)

    # Keep only pairs near the expected beam geometry
    dx_med = df["dx"].median()
    dy_med = df["dy"].median()
    sep_med = df["sep"].median()

    print(f"Median dx = {dx_med:.3f}")
    print(f"Median dy = {dy_med:.3f}")
    print(f"Median sep = {sep_med:.3f}")

    # Conservative cuts
    good = df[
        (np.abs(df["dx"] - dx_med) < 2.0) &
        (np.abs(df["dy"] - dy_med) < 2.0) &
        (np.abs(df["sep"] - sep_med) < 2.5)
    ].copy()

    out_csv = POL_DIR / "beam_pairs_filtered.csv"
    good.to_csv(out_csv, index=False)

    print(f"Original pairs: {len(df)}")
    print(f"Filtered good pairs: {len(good)}")
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    filter_pairs()