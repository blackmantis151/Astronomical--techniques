from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import numpy as np
from astropy.io import fits


BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data_raw"
MASTER_DIR = BASE_DIR / "output" / "masters"
CALIB_DIR = BASE_DIR / "output" / "calibrated"


def list_fits(folder: Path) -> List[Path]:
    patterns = ["*.fts", "*.fits", "*.fit", "*.FTS", "*.FITS", "*.FIT"]
    files: List[Path] = []
    for pattern in patterns:
        files.extend(folder.glob(pattern))
    return sorted(files)


def read_fits(path: Path) -> Tuple[np.ndarray, fits.Header]:
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(np.float64)
        header = hdul[0].header.copy()
    return data, header


def main() -> None:
    parser = argparse.ArgumentParser(description="Calibrate science images.")
    parser.add_argument("--science-night", required=True, help="Science night, e.g. 20200824")
    args = parser.parse_args()

    science_night = args.science_night
    light_dir = RAW_DIR / science_night / "light"
    master_dir = MASTER_DIR / science_night
    out_dir = CALIB_DIR / science_night
    out_dir.mkdir(parents=True, exist_ok=True)

    light_files = list_fits(light_dir)
    if not light_files:
        raise FileNotFoundError(f"No light files found in {light_dir}")

    master_bias, _ = read_fits(master_dir / "master_bias.fts")
    master_dark, _ = read_fits(master_dir / "master_dark.fts")
    master_flat, _ = read_fits(master_dir / "master_flat.fts")

    print(f"Science night: {science_night}")
    print(f"Light files found: {len(light_files)}")

    for i, file in enumerate(light_files, start=1):
        data, header = read_fits(file)
        calibrated = (data - master_bias - master_dark) / master_flat

        header["HISTORY"] = "Calibrated using (raw - master_bias - master_dark) / master_flat"

        out_name = file.stem + "_calib.fts"
        out_path = out_dir / out_name

        fits.writeto(out_path, calibrated, header=header, overwrite=True)
        print(f"[{i}/{len(light_files)}] Saved {out_path.name}")

    print("\nAll calibrated files saved in:")
    print(out_dir)


if __name__ == "__main__":
    main()