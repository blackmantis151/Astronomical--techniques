from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import ImageNormalize, PercentileInterval, SqrtStretch


BASE_DIR = Path(__file__).resolve().parent.parent
MASTER_DIR = BASE_DIR / "output" / "masters"
CALIB_DIR = BASE_DIR / "output" / "calibrated"
PLOT_DIR = BASE_DIR / "output" / "plots"


def list_fits(folder: Path) -> List[Path]:
    patterns = ["*.fts", "*.fits", "*.fit", "*.FTS", "*.FITS", "*.FIT"]
    files = []
    for pattern in patterns:
        files.extend(folder.glob(pattern))
    return sorted(files)


def read_fits(path: Path):
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(float)
    return data


def save_png(data, outpath, title=""):
    interval = PercentileInterval(99.5)
    vmin, vmax = interval.get_limits(data)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())

    plt.figure(figsize=(7, 7))
    plt.imshow(data, origin="lower", cmap="gray", norm=norm)
    plt.colorbar(label="Counts")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Make PNG previews of master and calibrated FITS files.")
    parser.add_argument("--science-night", required=True, help="Science night, e.g. 20200824")
    args = parser.parse_args()

    night = args.science_night

    master_night_dir = MASTER_DIR / night
    calib_night_dir = CALIB_DIR / night

    out_master_dir = PLOT_DIR / "masters" / night
    out_calib_dir = PLOT_DIR / "calibrated" / night
    out_master_dir.mkdir(parents=True, exist_ok=True)
    out_calib_dir.mkdir(parents=True, exist_ok=True)

    master_files = list_fits(master_night_dir)
    calib_files = list_fits(calib_night_dir)

    print(f"Making master previews for {night} ...")
    for f in master_files:
        data = read_fits(f)
        outpath = out_master_dir / f"{f.stem}.png"
        save_png(data, outpath, title=f.name)
        print(f"Saved {outpath}")

    print(f"\nMaking calibrated previews for {night} ...")
    for f in calib_files:
        data = read_fits(f)
        outpath = out_calib_dir / f"{f.stem}.png"
        save_png(data, outpath, title=f.name)
        print(f"Saved {outpath}")

    print("\nDone.")


if __name__ == "__main__":
    main()