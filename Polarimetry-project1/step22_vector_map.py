import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

from config import STACK_DIR, PLOT_DIR, POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data


def plot_vector_map(scale=80.0):
    ensure_directories(OUTPUT_SUBDIRS)

    data, _ = load_fits_data(STACK_DIR / "group_4_stack.fits")
    df = pd.read_csv(POL_DIR / "polarimetry_results_good.csv")

    norm = simple_norm(data, stretch="linear", percent=99.5)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(data, origin="lower", cmap="gray", norm=norm)

    for _, row in df.iterrows():
        x = row["x"]
        y = row["y"]
        P = row["P"]
        theta = np.radians(row["theta_deg"])

        L = scale * P
        dx = L * np.cos(theta)
        dy = L * np.sin(theta)

        ax.plot([x - dx/2, x + dx/2], [y - dy/2, y + dy/2], color="red", lw=1.5)
        ax.plot(x, y, marker="o", color="cyan", markersize=2)

    ax.set_xlabel("X pixel")
    ax.set_ylabel("Y pixel")
    ax.set_title("Polarization vector map")

    out_png = PLOT_DIR / "polarization_vector_map.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out_png}")
    print(f"Plotted {len(df)} good polarization vectors")


if __name__ == "__main__":
    plot_vector_map()