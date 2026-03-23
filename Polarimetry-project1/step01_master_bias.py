import numpy as np
from astropy.io import fits

from config import BIAS_DIR, MASTER_BIAS_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data, save_fits, show_and_save_image


def make_master_bias():
    ensure_directories(OUTPUT_SUBDIRS)

    bias_files = sorted(BIAS_DIR.glob("*.fits"))
    if len(bias_files) == 0:
        raise FileNotFoundError(f"No bias files found in {BIAS_DIR}")

    stack = []

    for f in bias_files:
        data, _ = load_fits_data(f)
        stack.append(data)

    stack = np.array(stack, dtype=float)
    master_bias = np.median(stack, axis=0)

    clean_header = fits.Header()
    clean_header["IMAGETYP"] = "MASTER_BIAS"
    clean_header["NCOMBINE"] = len(bias_files)
    clean_header["COMMENT"] = "Median combined master bias"

    out_fits = MASTER_BIAS_DIR / "master_bias.fits"
    out_png = MASTER_BIAS_DIR / "master_bias.png"

    save_fits(out_fits, master_bias, clean_header)
    show_and_save_image(master_bias, out_png, title="Master Bias")

    print(f"Saved master bias: {out_fits}")
    print(f"Saved plot: {out_png}")


if __name__ == "__main__":
    make_master_bias()