import pandas as pd
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

from config import STACK_DIR, PHOT_DIR, PLOT_DIR, DAO_FWHM, DAO_THRESHOLD_SIGMA, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data, plot_detected_sources


def detect_sources_in_reference_group(group_id=1):
    ensure_directories(OUTPUT_SUBDIRS)

    stack_path = STACK_DIR / f"group_{group_id}_stack.fits"
    data, _ = load_fits_data(stack_path)

    mean, med, std = sigma_clipped_stats(data, sigma=3.0)
    finder = DAOStarFinder(fwhm=DAO_FWHM, threshold=DAO_THRESHOLD_SIGMA * std)
    sources = finder(data - med)

    if sources is None:
        raise RuntimeError("No sources detected")

    src_df = sources.to_pandas()

    out_csv = PHOT_DIR / f"group_{group_id}_sources.csv"
    out_png = PLOT_DIR / f"group_{group_id}_detected_sources.png"

    src_df.to_csv(out_csv, index=False)
    plot_detected_sources(
        data,
        src_df,
        outpath=out_png,
        title=f"Detected sources in group {group_id} stack",
        radius=6,
    )

    print(f"Detected {len(src_df)} sources")
    print(f"Saved source table: {out_csv}")
    print(f"Saved source plot: {out_png}")


if __name__ == "__main__":
    detect_sources_in_reference_group()