import numpy as np
import pandas as pd

from config import PHOT_DIR, POL_DIR, PLOT_DIR, STACK_DIR, DX_BEAM, DY_BEAM, PAIR_TOL, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data, plot_beam_pairs


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
                "dx": xy[j, 0] - x,
                "dy": xy[j, 1] - y,
                "sep": np.sqrt((xy[j, 0] - x)**2 + (xy[j, 1] - y)**2),
            })
            used.add(i)
            used.add(j)

    pair_df = pd.DataFrame(pairs)

    out_csv = POL_DIR / "beam_pairs.csv"
    pair_df.to_csv(out_csv, index=False)

    data, _ = load_fits_data(STACK_DIR / "group_1_stack.fits")
    out_png = PLOT_DIR / "beam_pairs_overlay.png"

    plot_beam_pairs(
        data,
        pair_df,
        outpath=out_png,
        title="Beam pair overlay on group 1 stack",
    )

    print(f"Paired {len(pair_df)} O/E sources")
    print(f"Saved pairs: {out_csv}")
    print(f"Saved overlay: {out_png}")


if __name__ == "__main__":
    pair_beams()