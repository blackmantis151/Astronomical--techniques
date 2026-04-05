import pandas as pd
import matplotlib.pyplot as plt

from config import PHOT_DIR, GROUP_TO_HWP, PLOT_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def plot_modulation_curve(example_pair_id=2, chosen_aperture=5):
    ensure_directories(OUTPUT_SUBDIRS)

    df = pd.read_csv(PHOT_DIR / "aperture_photometry.csv")
    s = df[(df["pair_id"] == example_pair_id) & (df["r_ap"] == chosen_aperture)].copy()

    s["hwp"] = s["group_id"].map(GROUP_TO_HWP)
    s = s.sort_values("hwp")

    s["ratio"] = s["flux_o"] / s["flux_e"]
    s["diffnorm"] = (s["flux_o"] - s["flux_e"]) / (s["flux_o"] + s["flux_e"])

    # O and E fluxes
    plt.figure(figsize=(8, 5))
    plt.plot(s["hwp"], s["flux_o"], marker="o", label="O beam")
    plt.plot(s["hwp"], s["flux_e"], marker="s", label="E beam")
    plt.xlabel("HWP angle (deg)")
    plt.ylabel("Flux")
    plt.title(f"O/E flux vs HWP angle for pair {example_pair_id}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    out1 = PLOT_DIR / f"modulation_flux_pair_{example_pair_id}.png"
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()

    # Ratio modulation
    plt.figure(figsize=(8, 5))
    plt.plot(s["hwp"], s["ratio"], marker="o")
    plt.xlabel("HWP angle (deg)")
    plt.ylabel("O/E flux ratio")
    plt.title(f"O/E ratio modulation for pair {example_pair_id}")
    plt.grid(True)
    plt.tight_layout()
    out2 = PLOT_DIR / f"modulation_ratio_pair_{example_pair_id}.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()

    # Normalized difference modulation
    plt.figure(figsize=(8, 5))
    plt.plot(s["hwp"], s["diffnorm"], marker="o")
    plt.xlabel("HWP angle (deg)")
    plt.ylabel("(O-E)/(O+E)")
    plt.title(f"Normalized modulation for pair {example_pair_id}")
    plt.grid(True)
    plt.tight_layout()
    out3 = PLOT_DIR / f"modulation_diffnorm_pair_{example_pair_id}.png"
    plt.savefig(out3, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out1}")
    print(f"Saved: {out2}")
    print(f"Saved: {out3}")


if __name__ == "__main__":
    plot_modulation_curve(example_pair_id=2, chosen_aperture=5)