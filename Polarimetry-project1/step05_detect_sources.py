import pandas as pd
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

from config import STACK_DIR, PHOT_DIR, DAO_FWHM, DAO_THRESHOLD_SIGMA, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data


def detect_sources_in_reference_group(group_id=1):
    ensure_directories(OUTPUT_SUBDIRS)

    stack_path = STACK_DIR / f"group_{group_id}_stack.fits"
    data, _ = load_fits_data(stack_path)

    mean, med, std = sigma_clipped_stats(data, sigma=3.0)
    finder = DAOStarFinder(fwhm=DAO_FWHM, threshold=DAO_THRESHOLD_SIGMA * std)
    sources = finder(data - med)

    if sources is None:
        raise RuntimeError("No sources detected")

    out_csv = PHOT_DIR / f"group_{group_id}_sources.csv"
    sources.to_pandas().to_csv(out_csv, index=False)

    print(f"Detected {len(sources)} sources")
    print(f"Saved source table: {out_csv}")


if __name__ == "__main__":
    detect_sources_in_reference_group()