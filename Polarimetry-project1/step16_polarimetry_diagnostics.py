import pandas as pd
import matplotlib.pyplot as plt

from config import POL_DIR, PLOT_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def make_polarimetry_diagnostics():
    ensure_directories(OUTPUT_SUBDIRS)

    df = pd.read_csv(POL_DIR / "polarimetry_results.csv")

    # q-u plane
    plt.figure(figsize=(6, 6))
    plt.scatter(df["q"], df["u"])
    plt.axhline(0, color="black", lw=1)
    plt.axvline(0, color="black", lw=1)
    plt.xlabel("q")
    plt.ylabel("u")
    plt.title("q-u plane")
    plt.grid(True)
    plt.tight_layout()
    out1 = PLOT_DIR / "qu_plane.png"
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()

    # P histogram
    plt.figure(figsize=(7, 5))
    plt.hist(100 * df["P"], bins=15)
    plt.xlabel("Polarization (%)")
    plt.ylabel("Number of stars")
    plt.title("Polarization fraction distribution")
    plt.tight_layout()
    out2 = PLOT_DIR / "polarization_histogram.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()

    # Theta histogram
    plt.figure(figsize=(7, 5))
    plt.hist(df["theta_deg"], bins=15)
    plt.xlabel("Polarization angle (deg)")
    plt.ylabel("Number of stars")
    plt.title("Polarization angle distribution")
    plt.tight_layout()
    out3 = PLOT_DIR / "polarization_angle_histogram.png"
    plt.savefig(out3, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out1}")
    print(f"Saved: {out2}")
    print(f"Saved: {out3}")


if __name__ == "__main__":
    make_polarimetry_diagnostics()