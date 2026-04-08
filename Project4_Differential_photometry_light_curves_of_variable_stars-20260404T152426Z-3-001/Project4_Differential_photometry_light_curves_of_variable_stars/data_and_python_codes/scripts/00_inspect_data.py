#!/usr/bin/env python3

from pathlib import Path
from astropy.io import fits
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data_raw"
OUTPUT_DIR = BASE_DIR / "output" / "logs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

KEYWORDS = [
    "OBJECT",
    "FILTER",
    "EXPTIME",
    "IMAGETYP",
    "DATE-OBS",
    "JD",
    "MJD",
    "NAXIS1",
    "NAXIS2",
]

def list_fits_files(folder: Path):
    patterns = ["*.fts", "*.fits", "*.fit", "*.FTS", "*.FITS", "*.FIT"]
    files = []
    for pattern in patterns:
        files.extend(folder.rglob(pattern))
    return sorted(files)

def safe_get(header, key):
    return header[key] if key in header else None

def classify_from_path(filepath: Path):
    parts = [p.lower() for p in filepath.parts]
    if "bias" in parts:
        return "bias"
    if "dark" in parts:
        return "dark"
    if "flat" in parts:
        return "flat"
    if "light" in parts:
        return "light"
    return "unknown"

def get_night_from_path(filepath: Path):
    for part in filepath.parts:
        if part.isdigit() and len(part) == 8:
            return part
    return None

def inspect_fits_file(filepath: Path):
    try:
        with fits.open(filepath) as hdul:
            header = hdul[0].header
            data = hdul[0].data

            row = {
                "file": str(filepath.relative_to(BASE_DIR)),
                "night": get_night_from_path(filepath),
                "type_from_path": classify_from_path(filepath),
                "shape": None if data is None else str(data.shape),
            }

            for key in KEYWORDS:
                row[key] = safe_get(header, key)

            return row

    except Exception as e:
        return {
            "file": str(filepath.relative_to(BASE_DIR)),
            "night": get_night_from_path(filepath),
            "type_from_path": classify_from_path(filepath),
            "error": str(e),
        }

def main():
    print(f"BASE_DIR   : {BASE_DIR}")
    print(f"RAW_DIR    : {RAW_DIR}")
    print(f"RAW exists : {RAW_DIR.exists()}")

    fits_files = list_fits_files(RAW_DIR)

    print(f"Found files: {len(fits_files)}")

    if not fits_files:
        print("No FITS-like files found.")
        return

    rows = [inspect_fits_file(f) for f in fits_files]
    df = pd.DataFrame(rows)

    csv_path = OUTPUT_DIR / "fits_inventory.csv"
    df.to_csv(csv_path, index=False)

    print("\nSaved FITS inventory to:")
    print(csv_path)

    print("\nCounts by night and type:")
    print(df.groupby(["night", "type_from_path"]).size())

    if "FILTER" in df.columns:
        filters = df["FILTER"].dropna().unique()
        print("\nUnique filters:")
        print(filters)

    if "EXPTIME" in df.columns:
        exptimes = sorted(df["EXPTIME"].dropna().unique())
        print("\nUnique exposure times:")
        print(exptimes)

    print("\nFirst 20 rows:")
    print(df.head(20).to_string(index=False))

if __name__ == "__main__":
    main()