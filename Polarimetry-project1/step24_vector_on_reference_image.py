import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.visualization import simple_norm

from config import STACK_DIR, POL_DIR, PLOT_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data


# =========================================
# SETTINGS
# =========================================

# Use the observed stacked image as background
OBS_IMAGE_PATH = STACK_DIR / "group_1_stack.fits"

# Use good-quality results if available
POL_TABLE = POL_DIR / "polarimetry_results.csv"
# If you have a cleaned table, use this instead:
# POL_TABLE = POL_DIR / "polarimetry_results_good.csv"

OUT_PNG = PLOT_DIR / "polarization_vectors_on_oriented_observed.png"

VECTOR_SCALE = 120.0
MIN_P = 0.0


# =========================================
# IMAGE ORIENTATION
# =========================================

def orient_image(data):
    """
    Apply the same orientation correction recommended for matching sky images:
    1. invert X
    2. rotate 90 degrees anticlockwise
    """
    data = np.fliplr(data)
    data = np.rot90(data, k=1)
    return data


def transform_xy_for_oriented_image(x, y, nx, ny):
    """
    Transform coordinates from original detector frame
    into the oriented image frame using the same operation:
    1. invert X
    2. rotate 90 degrees anticlockwise
    """
    # step 1: invert X
    x1 = (nx - 1) - x
    y1 = y

    # step 2: rotate 90 deg anticlockwise
    x2 = y1
    y2 = (nx - 1) - x1

    return x2, y2


# =========================================
# MAIN
# =========================================

def main():
    ensure_directories(OUTPUT_SUBDIRS)

    # Load original observed image
    obs_data, _ = load_fits_data(OBS_IMAGE_PATH)
    ny, nx = obs_data.shape

    # Create oriented version for plotting
    oriented_data = orient_image(obs_data)

    # Load polarization table
    df = pd.read_csv(POL_TABLE).copy()

    # Keep only finite valid rows
    df = df[
        np.isfinite(df["x"]) &
        np.isfinite(df["y"]) &
        np.isfinite(df["P"]) &
        np.isfinite(df["theta_deg"])
    ].copy()

    df = df[df["P"] >= MIN_P].copy()

    # Transform coordinates into oriented image frame
    x_plot = []
    y_plot = []

    for _, row in df.iterrows():
        xt, yt = transform_xy_for_oriented_image(row["x"], row["y"], nx, ny)
        x_plot.append(xt)
        y_plot.append(yt)

    df["x_plot"] = x_plot
    df["y_plot"] = y_plot

    # Plot background
    norm = simple_norm(oriented_data, stretch="linear", percent=99.5)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(oriented_data, origin="lower", cmap="gray", norm=norm)

    # Plot vectors
    for _, row in df.iterrows():
        x = row["x_plot"]
        y = row["y_plot"]
        P = row["P"]
        theta = np.radians(row["theta_deg"])

        # Polarization is headless, so line segment centered on source is fine
        L = VECTOR_SCALE * P
        dx = L * np.cos(theta)
        dy = L * np.sin(theta)

        ax.plot(
            [x - dx / 2, x + dx / 2],
            [y - dy / 2, y + dy / 2],
            color="red",
            lw=1.5
        )

        ax.plot(x, y, marker="o", color="cyan", markersize=2)

    ax.set_xlabel("X pixel")
    ax.set_ylabel("Y pixel")
    ax.set_title("Polarization vectors on oriented observed image")

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {OUT_PNG}")
    print(f"Plotted {len(df)} vectors")


if __name__ == "__main__":
    main()