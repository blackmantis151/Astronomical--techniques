import pandas as pd
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

from config import POL_DIR, STACK_DIR, PLOT_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data


def plot_filtered_pairs_on_image(data, pair_df, outpath, title="Filtered beam pairs", cmap="gray"):
    """
    Plot filtered O/E beam pairs on the stacked image.

    Yellow = O beam
    Red = E beam
    Cyan line = connection between the pair
    """
    norm = simple_norm(data, stretch="linear", percent=99.5)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(data, origin="lower", cmap=cmap, norm=norm)

    for _, row in pair_df.iterrows():
        x_o = row["x_o"]
        y_o = row["y_o"]
        x_e = row["x_e"]
        y_e = row["y_e"]

        ax.plot([x_o, x_e], [y_o, y_e], color="cyan", linewidth=1.2)
        ax.plot(x_o, y_o, marker="o", color="yellow", markersize=4)
        ax.plot(x_e, y_e, marker="o", color="red", markersize=4)

    ax.set_xlabel("X pixel")
    ax.set_ylabel("Y pixel")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()


def main_step18():
    ensure_directories(OUTPUT_SUBDIRS)

    pair_csv = POL_DIR / "beam_pairs_filtered.csv"
    stack_fits = STACK_DIR / "group_1_stack.fits"
    out_png = PLOT_DIR / "beam_pairs_filtered_overlay.png"

    pair_df = pd.read_csv(pair_csv)
    data, _ = load_fits_data(stack_fits)

    plot_filtered_pairs_on_image(
        data=data,
        pair_df=pair_df,
        outpath=out_png,
        title="Filtered beam pairs on group 1 stack",
    )

    print(f"Loaded filtered pairs: {pair_csv}")
    print(f"Number of filtered pairs: {len(pair_df)}")
    print(f"Saved overlay: {out_png}")


if __name__ == "__main__":
    main_step18()