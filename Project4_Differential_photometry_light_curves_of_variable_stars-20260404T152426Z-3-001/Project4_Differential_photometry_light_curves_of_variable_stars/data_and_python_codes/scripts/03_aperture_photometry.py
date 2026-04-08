

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS


BASE_DIR = Path(__file__).resolve().parent.parent
CALIB_DIR = BASE_DIR / "output" / "calibrated"
PHOT_DIR = BASE_DIR / "output" / "photometry"
PLOT_DIR = BASE_DIR / "output" / "plots"
REGION_FILE = BASE_DIR / "RefStar_BLLac.reg"

PHOT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

REFERENCE_MAGS = {
    "BL-Lac": np.nan,
    "ref-A": 13.73,
    "ref-B": 11.99,
    "ref-C": 13.79,
}


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


def parse_ds9_region_file(path: Path) -> Dict[str, Dict]:
    if not path.exists():
        raise FileNotFoundError(f"Region file not found: {path}")

    circle_re = re.compile(
        r"circle\(([^,]+),([^,]+),([^)]+)\).*text=\{([^}]+)\}",
        re.IGNORECASE,
    )
    annulus_re = re.compile(
        r"annulus\(([^,]+),([^,]+),([^,]+),([^)]+)\).*text=\{([^}]+)\}",
        re.IGNORECASE,
    )

    circles = {}
    annuli = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            c_match = circle_re.search(line)
            if c_match:
                ra_str, dec_str, r_str, label = c_match.groups()
                circles[label] = {
                    "ra_str": ra_str.strip(),
                    "dec_str": dec_str.strip(),
                    "r_arcsec": float(r_str.replace('"', "").strip()),
                }
                continue

            a_match = annulus_re.search(line)
            if a_match:
                ra_str, dec_str, rin_str, rout_str, label = a_match.groups()
                annuli[label] = {
                    "ra_str": ra_str.strip(),
                    "dec_str": dec_str.strip(),
                    "rin_arcsec": float(rin_str.replace('"', "").strip()),
                    "rout_arcsec": float(rout_str.replace('"', "").strip()),
                }

    labels = sorted(set(circles) & set(annuli))
    if not labels:
        raise ValueError("No matching circle and annulus entries found in the region file.")

    objects = {}
    for label in labels:
        coord = SkyCoord(
            circles[label]["ra_str"],
            circles[label]["dec_str"],
            unit=(u.hourangle, u.deg),
            frame="icrs",
        )
        objects[label] = {
            "coord": coord,
            "r_arcsec": circles[label]["r_arcsec"],
            "rin_arcsec": annuli[label]["rin_arcsec"],
            "rout_arcsec": annuli[label]["rout_arcsec"],
        }

    return objects


def circular_mask(shape: Tuple[int, int], x0: float, y0: float, r: float) -> np.ndarray:
    yy, xx = np.indices(shape)
    return (xx - x0) ** 2 + (yy - y0) ** 2 <= r ** 2


def annulus_mask(shape: Tuple[int, int], x0: float, y0: float, r_in: float, r_out: float) -> np.ndarray:
    yy, xx = np.indices(shape)
    rr2 = (xx - x0) ** 2 + (yy - y0) ** 2
    return (rr2 >= r_in ** 2) & (rr2 <= r_out ** 2)


def measure_flux(
    data: np.ndarray,
    x: float,
    y: float,
    r_pix: float,
    rin_pix: float,
    rout_pix: float,
) -> Tuple[float, float]:
    src_mask = circular_mask(data.shape, x, y, r_pix)
    sky_mask = annulus_mask(data.shape, x, y, rin_pix, rout_pix)

    source_pixels = data[src_mask]
    sky_pixels = data[sky_mask]

    if sky_pixels.size == 0:
        raise ValueError("Sky annulus has zero pixels.")
    if source_pixels.size == 0:
        raise ValueError("Source aperture has zero pixels.")

    _, sky_median, sky_std = sigma_clipped_stats(sky_pixels, sigma=3.0)
    net_flux = np.sum(source_pixels) - sky_median * source_pixels.size

    return float(net_flux), float(sky_std)


def get_pixel_scale_arcsec(header, wcs, file):
    try:
        pixel_scales = wcs.proj_plane_pixel_scales()
        pixel_scale_arcsec = np.mean([
            s.to(u.arcsec).value if hasattr(s, "to") else float(s) * 3600.0
            for s in pixel_scales
        ])
        return pixel_scale_arcsec
    except Exception:
        cdelt1 = header.get("CDELT1")
        cdelt2 = header.get("CDELT2")
        if cdelt1 is not None and cdelt2 is not None:
            return np.mean([abs(cdelt1), abs(cdelt2)]) * 3600.0
        raise ValueError(f"Could not determine WCS pixel scale for {file}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Perform differential aperture photometry.")
    parser.add_argument(
        "--nights",
        nargs="+",
        required=True,
        help="Science nights, e.g. 20200824 20200825",
    )
    args = parser.parse_args()

    region_objects = parse_ds9_region_file(REGION_FILE)
    all_rows = []

    for night in args.nights:
        night_dir = CALIB_DIR / night
        files = list_fits(night_dir)
        if not files:
            raise FileNotFoundError(f"No calibrated files found in {night_dir}")

        print(f"\nProcessing night {night}")
        print(f"Number of calibrated files: {len(files)}")

        for file in files:
            data, header = read_fits(file)
            wcs = WCS(header)

            jd = header.get("JD")
            date_obs = header.get("DATE-OBS")
            obj_name = header.get("OBJECT", "UNKNOWN")

            pixel_scale_arcsec = get_pixel_scale_arcsec(header, wcs, file)

            row = {
                "night": night,
                "file": file.name,
                "JD": jd,
                "DATE-OBS": date_obs,
                "OBJECT": obj_name,
            }

            inst_mags = {}
            zero_points = []

            for label, info in region_objects.items():
                x_pix, y_pix = wcs.world_to_pixel(info["coord"])
                x_pix = float(x_pix)
                y_pix = float(y_pix)

                r_pix = info["r_arcsec"] / pixel_scale_arcsec
                rin_pix = info["rin_arcsec"] / pixel_scale_arcsec
                rout_pix = info["rout_arcsec"] / pixel_scale_arcsec

                flux, sky_std = measure_flux(data, x_pix, y_pix, r_pix, rin_pix, rout_pix)

                row[f"{label}_xpix"] = x_pix
                row[f"{label}_ypix"] = y_pix
                row[f"{label}_flux"] = flux
                row[f"{label}_sky_std"] = sky_std

                inst_mag = -2.5 * np.log10(flux) if flux > 0 else np.nan
                row[f"{label}_instmag"] = inst_mag
                inst_mags[label] = inst_mag

                ref_mag = REFERENCE_MAGS.get(label, np.nan)
                if np.isfinite(ref_mag) and np.isfinite(inst_mag):
                    zero_points.append(ref_mag - inst_mag)

            if zero_points and np.isfinite(inst_mags.get("BL-Lac", np.nan)):
                zero_point = float(np.median(zero_points))
                row["zero_point"] = zero_point
                row["target_Rmag"] = inst_mags["BL-Lac"] + zero_point
            else:
                row["zero_point"] = np.nan
                row["target_Rmag"] = np.nan

            ref_inst_mags = [
                inst_mags.get("ref-A", np.nan),
                inst_mags.get("ref-B", np.nan),
                inst_mags.get("ref-C", np.nan),
            ]
            ref_inst_mags = [m for m in ref_inst_mags if np.isfinite(m)]

            if np.isfinite(inst_mags.get("BL-Lac", np.nan)) and ref_inst_mags:
                row["target_minus_ref_ensemble"] = inst_mags["BL-Lac"] - float(np.mean(ref_inst_mags))
            else:
                row["target_minus_ref_ensemble"] = np.nan

            all_rows.append(row)

    df = pd.DataFrame(all_rows)
    if "JD" in df.columns:
        df = df.sort_values(by="JD", na_position="last").reset_index(drop=True)

    night_label = "-".join(args.nights)
    csv_path = PHOT_DIR / f"BLLac_{night_label}.csv"
    df.to_csv(csv_path, index=False)

    plt.figure(figsize=(8, 5))
    plt.plot(df["JD"], df["target_Rmag"], "o")
    plt.gca().invert_yaxis()
    plt.xlabel("JD")
    plt.ylabel("Target R magnitude")
    plt.title(f"BL Lac Light Curve ({night_label})")
    plt.tight_layout()
    mag_plot = PLOT_DIR / f"BLLac_Rmag_{night_label}.png"
    plt.savefig(mag_plot, dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["JD"], df["target_minus_ref_ensemble"], "o")
    plt.xlabel("JD")
    plt.ylabel("Differential instrumental magnitude")
    plt.title(f"BL Lac Differential Light Curve ({night_label})")
    plt.tight_layout()
    diff_plot = PLOT_DIR / f"BLLac_diffmag_{night_label}.png"
    plt.savefig(diff_plot, dpi=200)
    plt.close()

    print("\nSaved photometry table:")
    print(csv_path)
    print("\nSaved plots:")
    print(mag_plot)
    print(diff_plot)


if __name__ == "__main__":
    main()