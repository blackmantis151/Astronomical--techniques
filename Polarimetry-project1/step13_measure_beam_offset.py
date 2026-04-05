import numpy as np
import matplotlib.pyplot as plt

from astropy.visualization import simple_norm

from config import STACK_DIR, PLOT_DIR, POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data


def measure_beam_offset():
    ensure_directories(OUTPUT_SUBDIRS)

    data, _ = load_fits_data(STACK_DIR / "group_1_stack.fits")

    # Edit these after inspecting the central double source
    # Current trial values
    x_o, y_o = 205, 216
    x_e, y_e = 236, 221

    dx = x_e - x_o
    dy = y_e - y_o
    sep = np.sqrt(dx**2 + dy**2)

    # Save values
    out_txt = POL_DIR / "beam_offset.txt"
    with open(out_txt, "w") as f:
        f.write("Beam offset from central double source\n")
        f.write("======================================\n")
        f.write(f"x_o = {x_o}\n")
        f.write(f"y_o = {y_o}\n")
        f.write(f"x_e = {x_e}\n")
        f.write(f"y_e = {y_e}\n")
        f.write(f"dx = {dx}\n")
        f.write(f"dy = {dy}\n")
        f.write(f"sep = {sep}\n")

    # Full image plot
    norm = simple_norm(data, stretch="linear", percent=99.5)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(data, origin="lower", cmap="gray", norm=norm)
    ax.plot([x_o, x_e], [y_o, y_e], color="red", lw=2)
    ax.plot(x_o, y_o, "yo", markersize=7, label="O beam")
    ax.plot(x_e, y_e, "co", markersize=7, label="E beam")
    ax.set_xlabel("X pixel")
    ax.set_ylabel("Y pixel")
    ax.set_title("Central double source beam offset")
    ax.legend()
    out_full = PLOT_DIR / "beam_offset_central_source.png"
    plt.savefig(out_full, dpi=150, bbox_inches="tight")
    plt.close()

    # Zoomed cutout
    cut_half = 25
    x1 = max(0, int(min(x_o, x_e) - cut_half))
    x2 = min(data.shape[1], int(max(x_o, x_e) + cut_half))
    y1 = max(0, int(min(y_o, y_e) - cut_half))
    y2 = min(data.shape[0], int(max(y_o, y_e) + cut_half))

    sub = data[y1:y2, x1:x2]
    norm_sub = simple_norm(sub, stretch="linear", percent=99.5)

    fig2, ax2 = plt.subplots(figsize=(6, 6))
    ax2.imshow(sub, origin="lower", cmap="gray", norm=norm_sub)
    ax2.plot([x_o - x1, x_e - x1], [y_o - y1, y_e - y1], color="red", lw=2)
    ax2.plot(x_o - x1, y_o - y1, "yo", markersize=7, label="O beam")
    ax2.plot(x_e - x1, y_e - y1, "co", markersize=7, label="E beam")
    ax2.set_xlabel("X pixel in cutout")
    ax2.set_ylabel("Y pixel in cutout")
    ax2.set_title("Zoomed central beam pair")
    ax2.legend()
    out_zoom = PLOT_DIR / "beam_offset_central_zoom.png"
    plt.savefig(out_zoom, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out_txt}")
    print(f"Saved: {out_full}")
    print(f"Saved: {out_zoom}")
    print(f"dx = {dx}, dy = {dy}, separation = {sep:.2f} pixels")


if __name__ == "__main__":
    measure_beam_offset()