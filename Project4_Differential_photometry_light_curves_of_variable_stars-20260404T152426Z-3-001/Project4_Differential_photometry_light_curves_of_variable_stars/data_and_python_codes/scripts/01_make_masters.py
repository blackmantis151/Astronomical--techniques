#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import numpy as np
from astropy.io import fits


BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data_raw"
OUT_DIR = BASE_DIR / "output" / "masters"


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


def mean_combine(files: List[Path]) -> Tuple[np.ndarray, fits.Header]:
    if not files:
        raise FileNotFoundError("No FITS files found.")
    stack = []
    first_header = None
    for file in files:
        data, header = read_fits(file)
        if first_header is None:
            first_header = header
        stack.append(data)
    return np.mean(np.stack(stack, axis=0), axis=0), first_header


def main() -> None:
    parser = argparse.ArgumentParser(description="Create master bias, dark and flat.")
    parser.add_argument("--science-night", required=True, help="Science night, e.g. 20200824")
    parser.add_argument(
        "--flat-nights",
        nargs="+",
        required=True,
        help="Flat nights, e.g. 20200819 20200830",
    )
    parser.add_argument(
        "--flat-filter",
        default="R",
        help="Filter keyword to select flats, default: R",
    )
    args = parser.parse_args()

    science_night = args.science_night
    flat_nights = args.flat_nights
    flat_filter = args.flat_filter.lower()

    science_dir = RAW_DIR / science_night
    bias_dir = science_dir / "bias"
    dark_dir = science_dir / "dark"

    out_dir = OUT_DIR / science_night
    out_dir.mkdir(parents=True, exist_ok=True)

    bias_files = list_fits(bias_dir)
    dark_files = list_fits(dark_dir)

    if not bias_files:
        raise FileNotFoundError(f"No bias files found in {bias_dir}")
    if not dark_files:
        raise FileNotFoundError(f"No dark files found in {dark_dir}")

    print(f"Science night: {science_night}")
    print(f"Bias files found : {len(bias_files)}")
    print(f"Dark files found : {len(dark_files)}")

    master_bias, bias_header = mean_combine(bias_files)

    dark_stack = []
    dark_header = None
    for file in dark_files:
        dark_data, dark_hdr = read_fits(file)
        if dark_header is None:
            dark_header = dark_hdr
        dark_stack.append(dark_data - master_bias)
    master_dark = np.mean(np.stack(dark_stack, axis=0), axis=0)

    flat_files_all: List[Path] = []
    for night in flat_nights:
        flat_dir = RAW_DIR / night / "flat"
        flat_files_all.extend(list_fits(flat_dir))

    if not flat_files_all:
        raise FileNotFoundError("No flat files found in the supplied flat nights.")

    flat_files = [
        f for f in flat_files_all
        if flat_filter in f.name.lower()
    ]

    if not flat_files:
        print("No flat files matched filter keyword in filename.")
        print("Using all flat files provided.")
        flat_files = flat_files_all

    print(f"Flat files selected: {len(flat_files)}")

    flat_stack = []
    flat_header = None
    for file in flat_files:
        flat_data, flat_hdr = read_fits(file)
        if flat_header is None:
            flat_header = flat_hdr
        corrected = flat_data - master_bias - master_dark
        flat_stack.append(corrected)

    mean_flat = np.mean(np.stack(flat_stack, axis=0), axis=0)
    flat_norm = np.mean(mean_flat)

    if not np.isfinite(flat_norm) or flat_norm == 0:
        raise ValueError("Flat normalization failed.")

    master_flat = mean_flat / flat_norm

    bias_header["HISTORY"] = "Master bias created by 01_make_masters.py"
    dark_header["HISTORY"] = "Master dark = mean(dark - master_bias)"
    flat_header["HISTORY"] = "Master flat = mean(flat - master_bias - master_dark) / mean"

    fits.writeto(out_dir / "master_bias.fts", master_bias, header=bias_header, overwrite=True)
    fits.writeto(out_dir / "master_dark.fts", master_dark, header=dark_header, overwrite=True)
    fits.writeto(out_dir / "master_flat.fts", master_flat, header=flat_header, overwrite=True)

    print("\nSaved master files:")
    print(out_dir / "master_bias.fts")
    print(out_dir / "master_dark.fts")
    print(out_dir / "master_flat.fts")


if __name__ == "__main__":
    main()