import numpy as np
import pandas as pd

from config import PHOT_DIR, POL_DIR, DX_BEAM, DY_BEAM, PAIR_TOL, OUTPUT_SUBDIRS
from utils import ensure_directories


def pair_beams():
    ensure_directories(OUTPUT_SUBDIRS)

    src_csv = PHOT_DIR / "group_1_sources.csv"
    src = pd.read_csv(src_csv)

    xy = src[["xcentroid", "ycentroid"]].to_numpy()
    used = set()
    pairs = []

    for i, (x, y) in enumerate(xy):
        if i in used:
            continue

        target = np.array([x + DX_BEAM, y + DY_BEAM])
        dist = np.sqrt(np.sum((xy - target) ** 2, axis=1))
        j = np.argmin(dist)

        if dist[j] < PAIR_TOL and j != i and j not in used:
            pairs.append({
                "pair_id": len(pairs),
                "x_o": x,
                "y_o": y,
                "x_e": xy[j, 0],
                "y_e": xy[j, 1],
            })
            used.add(i)
            used.add(j)

    pair_df = pd.DataFrame(pairs)
    out_csv = POL_DIR / "beam_pairs.csv"
    pair_df.to_csv(out_csv, index=False)

    print(f"Paired {len(pair_df)} O/E sources")
    print(f"Saved pairs: {out_csv}")


if __name__ == "__main__":
    pair_beams()