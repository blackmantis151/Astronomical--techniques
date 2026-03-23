import numpy as np
from astropy.io import fits

from config import BIAS_SUB_DIR, STACK_DIR, PLOT_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data, save_fits, show_and_save_image, build_science_dataframe


def stack_groups():
    ensure_directories(OUTPUT_SUBDIRS)

    corrected_files = sorted(BIAS_SUB_DIR.glob("*.fits"))
    df = build_science_dataframe(corrected_files)

    for g in sorted(df["group_id"].dropna().unique()):
        group_files = df[df["group_id"] == g]["file"].tolist()
        stack = []

        for fname in group_files:
            data, _ = load_fits_data(BIAS_SUB_DIR / fname)
            stack.append(data)

        stack = np.array(stack, dtype=float)
        combined = np.median(stack, axis=0)

        clean_header = fits.Header()
        clean_header["IMAGETYP"] = "STACK"
        clean_header["GROUPID"] = int(g)
        clean_header["NCOMBINE"] = len(group_files)
        clean_header["COMMENT"] = "Median stack of bias-subtracted science frames"

        out_fits = STACK_DIR / f"group_{g}_stack.fits"
        out_png = PLOT_DIR / f"group_{g}_stack.png"

        save_fits(out_fits, combined, clean_header)
        show_and_save_image(combined, out_png, title=f"Group {g} Stack")

        print(f"Saved {out_fits}")


if __name__ == "__main__":
    stack_groups()