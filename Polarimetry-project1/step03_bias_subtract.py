from config import SCIENCE_DIR, MASTER_BIAS_DIR, BIAS_SUB_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data, save_fits, build_science_dataframe


def bias_subtract_science():
    ensure_directories(OUTPUT_SUBDIRS)

    master_bias_path = MASTER_BIAS_DIR / "master_bias.fits"
    master_bias, _ = load_fits_data(master_bias_path)

    science_files = sorted(SCIENCE_DIR.glob("*.fits"))
    df = build_science_dataframe(science_files)

    for _, row in df.iterrows():
        inpath = SCIENCE_DIR / row["file"]
        data, header = load_fits_data(inpath)

        corrected = data - master_bias
        outpath = BIAS_SUB_DIR / row["file"]
        save_fits(outpath, corrected, header)

    print(f"Bias-subtracted {len(df)} science files into {BIAS_SUB_DIR}")


if __name__ == "__main__":
    bias_subtract_science()