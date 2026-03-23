import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from config import (
    STACK_DIR,
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
    plot_marked_star,
    plot_cutout_with_circle,
    show_and_save_image
)


def mark_alignment_star():
    ensure_directories(OUTPUT_SUBDIRS)

    stack_path = STACK_DIR / f"group_{ALIGNMENT_GROUP_ID}_stack.fits"
    data, _ = load_fits_data(stack_path)

    full_out = PLOT_DIR / "alignment_star_full.png"
    cutout_out = PLOT_DIR / "alignment_star_cutout.png"

    plot_marked_star(
        data,
        ALIGNMENT_X,
        ALIGNMENT_Y,
        radius=ALIGNMENT_RADIUS,
        title=f"Alignment star on group {ALIGNMENT_GROUP_ID} stack",
        outpath=full_out,
    )

    plot_cutout_with_circle(
        data,
        ALIGNMENT_X,
        ALIGNMENT_Y,
        radius=ALIGNMENT_RADIUS,
        half_size=ALIGNMENT_CUTOUT_HALF_SIZE,
        title="Zoomed cutout of alignment star",
        outpath=cutout_out,
    )

    print(f"Saved full image check: {full_out}")
    print(f"Saved cutout check: {cutout_out}")


if __name__ == "__main__":
    mark_alignment_star()

