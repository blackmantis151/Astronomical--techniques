import pandas as pd
import matplotlib.pyplot as plt

from config import PHOT_DIR, PLOT_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def plot_photometry_diagnostics(example_pair_id=0):
    ensure_directories(OUTPUT_SUBDIRS)

    df = pd.read_csv(PHOT_DIR / "aperture_photometry.csv")
    s = df[df["pair_id"] == example_pair_id].copy()

    # Flux vs aperture
    plt.figure(figsize=(8, 5))
    for g in sorted(s["group_id"].unique()):
        sg = s[s["group_id"] == g].sort_values("r_ap")
        plt.plot(sg["r_ap"], sg["flux_o"], marker="o", label=f"O, group {g}")
        plt.plot(sg["r_ap"], sg["flux_e"], marker="s", linestyle="--", label=f"E, group {g}")

    plt.xlabel("Aperture radius (pixels)")
    plt.ylabel("Sky-subtracted flux")
    plt.title(f"Flux vs aperture for pair {example_pair_id}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    out1 = PLOT_DIR / f"photometry_flux_vs_aperture_pair_{example_pair_id}.png"
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()

    # Sky vs aperture
    plt.figure(figsize=(8, 5))
    for g in sorted(s["group_id"].unique()):
        sg = s[s["group_id"] == g].sort_values("r_ap")
        plt.plot(sg["r_ap"], sg["sky_o"], marker="o", label=f"Sky O, group {g}")
        plt.plot(sg["r_ap"], sg["sky_e"], marker="s", linestyle="--", label=f"Sky E, group {g}")

    plt.xlabel("Aperture radius (pixels)")
    plt.ylabel("Median sky per pixel")
    plt.title(f"Sky vs aperture for pair {example_pair_id}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    out2 = PLOT_DIR / f"photometry_sky_vs_aperture_pair_{example_pair_id}.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out1}")
    print(f"Saved: {out2}")


if __name__ == "__main__":
    plot_photometry_diagnostics(example_pair_id=0)