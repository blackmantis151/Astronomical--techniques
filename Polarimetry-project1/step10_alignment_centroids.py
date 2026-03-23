import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from photutils.centroids import centroid_2dg
from astropy.visualization import simple_norm

from config import (
    SCIENCE_DIR,
    PLOT_DIR,
    OUTPUT_SUBDIRS,
    ALIGNMENT_GROUP_ID,
    ALIGNMENT_X,
    ALIGNMENT_Y,
    ALIGNMENT_RADIUS,
    ALIGNMENT_CUTOUT_HALF_SIZE,
)
from utils import (
    ensure_directories,
    load_fits_data,
    build_science_dataframe,
    extract_cutout,
)


def measure_centroid_on_cutout(data, x, y, half_size):
    cutout, x1, y1 = extract_cutout(data, x, y, half_size=half_size)
    yc, xc = centroid_2dg(cutout)

    x_global = x1 + xc
    y_global = y1 + yc

    return cutout, x1, y1, xc, yc, x_global, y_global


def run_alignment_centroid_check():
    ensure_directories(OUTPUT_SUBDIRS)

    science_files = sorted(SCIENCE_DIR.glob("*.fits"))
    df = build_science_dataframe(science_files)
    df = df[df["group_id"] == ALIGNMENT_GROUP_ID].copy()

    rows = []

    for _, row in df.iterrows():
        fname = row["file"]
        data, _ = load_fits_data(SCIENCE_DIR / fname)

        cutout, x1, y1, xc, yc, xg, yg = measure_centroid_on_cutout(
            data,
            ALIGNMENT_X,
            ALIGNMENT_Y,
            ALIGNMENT_CUTOUT_HALF_SIZE,
        )

        norm = simple_norm(cutout, stretch="linear", percent=99.5)

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.imshow(cutout, origin="lower", cmap="gray", norm=norm)

        circ = Circle(
            (ALIGNMENT_X - x1, ALIGNMENT_Y - y1),
            ALIGNMENT_RADIUS,
            edgecolor="red",
            facecolor="none",
            linewidth=2,
        )
        ax.add_patch(circ)

        ax.plot(xc, yc, marker="+", color="yellow", markersize=12, mew=2)
        ax.set_title(f"{fname}")
        ax.set_xlabel("X pixel in cutout")
        ax.set_ylabel("Y pixel in cutout")

        outpath = PLOT_DIR / f"centroid_check_{fname.replace('.fits', '.png')}"
        plt.savefig(outpath, dpi=150, bbox_inches="tight")
        plt.close()

        rows.append({
            "file": fname,
            "group_id": ALIGNMENT_GROUP_ID,
            "x_guess": ALIGNMENT_X,
            "y_guess": ALIGNMENT_Y,
            "x_centroid": xg,
            "y_centroid": yg,
        })

    out_csv = PLOT_DIR / f"alignment_centroids_group_{ALIGNMENT_GROUP_ID}.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)

    print(f"Saved centroid table: {out_csv}")
    print(f"Saved centroid plots in: {PLOT_DIR}")


if __name__ == "__main__":
    run_alignment_centroid_check()