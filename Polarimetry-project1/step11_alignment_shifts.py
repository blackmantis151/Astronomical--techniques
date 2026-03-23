import pandas as pd
import matplotlib.pyplot as plt

from config import PLOT_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def plot_alignment_shifts(group_id=1):
    ensure_directories(OUTPUT_SUBDIRS)

    csv_path = PLOT_DIR / f"alignment_centroids_group_{group_id}.csv"
    df = pd.read_csv(csv_path)

    x_ref = df.iloc[0]["x_centroid"]
    y_ref = df.iloc[0]["y_centroid"]

    df["dx"] = df["x_centroid"] - x_ref
    df["dy"] = df["y_centroid"] - y_ref

    out_csv = PLOT_DIR / f"alignment_shifts_group_{group_id}.csv"
    df.to_csv(out_csv, index=False)

    plt.figure(figsize=(8, 5))
    plt.plot(df["file"], df["dx"], marker="o", label="dx")
    plt.plot(df["file"], df["dy"], marker="s", label="dy")
    plt.axhline(0, color="black", linewidth=1)
    plt.xticks(rotation=45)
    plt.ylabel("Shift (pixels)")
    plt.title(f"Alignment shifts for group {group_id}")
    plt.legend()
    plt.tight_layout()

    out_png = PLOT_DIR / f"alignment_shifts_group_{group_id}.png"
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved shift table: {out_csv}")
    print(f"Saved shift plot: {out_png}")


if __name__ == "__main__":
    plot_alignment_shifts(group_id=1)