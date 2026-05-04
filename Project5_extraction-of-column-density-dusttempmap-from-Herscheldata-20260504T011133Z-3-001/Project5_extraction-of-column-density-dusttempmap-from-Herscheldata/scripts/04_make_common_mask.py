from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_DIR = PROJECT_ROOT / "data" / "processed" / "reprojected_to_500grid"
OUT_DIR = PROJECT_ROOT / "data" / "processed" / "common_masked"
PLOT_DIR = PROJECT_ROOT / "plots" / "common_masked"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


MAP_FILES = {
    160: INPUT_DIR / "herschel_160um_mjysr_conv500_regrid500.fits",
    250: INPUT_DIR / "herschel_250um_mjysr_conv500_regrid500.fits",
    350: INPUT_DIR / "herschel_350um_mjysr_conv500_regrid500.fits",
    500: INPUT_DIR / "herschel_500um_mjysr_conv500_regrid500.fits",
}


def read_map(path):
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(float)
        header = hdul[0].header.copy()
    return data, header


def robust_limits(image):
    finite = np.isfinite(image)
    positive = finite & (image > 0)

    if positive.sum() > 100:
        values = image[positive]
    else:
        values = image[finite]

    if values.size == 0:
        return 0, 1

    vmin = np.nanpercentile(values, 2)
    vmax = np.nanpercentile(values, 99)

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin = np.nanmin(values)
        vmax = np.nanmax(values)

    return vmin, vmax


def save_map_plot(image, header, wavelength):
    fig = plt.figure(figsize=(7, 6))

    try:
        wcs = WCS(header)
        ax = plt.subplot(projection=wcs)
        im = ax.imshow(image, origin="lower", cmap="inferno")
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    except Exception:
        ax = plt.subplot()
        im = ax.imshow(image, origin="lower", cmap="inferno")
        ax.set_xlabel("Pixel x")
        ax.set_ylabel("Pixel y")

    vmin, vmax = robust_limits(image)
    im.set_clim(vmin, vmax)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Surface brightness [MJy/sr]")

    ax.set_title(f"{wavelength} micron common valid region")

    out_png = PLOT_DIR / f"common_{wavelength}um.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"Saved plot: {out_png}")


def save_mask_plot(mask, header):
    fig = plt.figure(figsize=(7, 6))

    try:
        wcs = WCS(header)
        ax = plt.subplot(projection=wcs)
        im = ax.imshow(mask.astype(float), origin="lower", cmap="gray")
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    except Exception:
        ax = plt.subplot()
        im = ax.imshow(mask.astype(float), origin="lower", cmap="gray")
        ax.set_xlabel("Pixel x")
        ax.set_ylabel("Pixel y")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Valid pixel mask")

    ax.set_title("Common valid mask for all Herschel bands")

    out_png = PLOT_DIR / "common_mask.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"Saved mask plot: {out_png}")


def main():
    print("\nMaking common valid mask for all Herschel bands\n")

    maps = {}
    headers = {}

    for wavelength, path in MAP_FILES.items():
        print("=" * 80)
        print(f"Reading {wavelength} micron map")
        print(f"File: {path}")

        if not path.exists():
            raise FileNotFoundError(f"Missing file: {path}")

        data, header = read_map(path)

        maps[wavelength] = data
        headers[wavelength] = header

        print(f"Shape: {data.shape}")
        print(f"Finite pixels: {np.isfinite(data).sum()} / {data.size}")
        print(f"Positive finite pixels: {(np.isfinite(data) & (data > 0)).sum()} / {data.size}")
        print(f"Min: {np.nanmin(data):.6e}")
        print(f"Max: {np.nanmax(data):.6e}")
        print(f"Median: {np.nanmedian(data):.6e}")

    shapes = [maps[w].shape for w in maps]
    if len(set(shapes)) != 1:
        raise ValueError(f"Maps do not have the same shape: {shapes}")

    common_mask = np.ones_like(maps[160], dtype=bool)

    for wavelength, data in maps.items():
        band_mask = np.isfinite(data) & (data > 0)
        common_mask &= band_mask

    print("=" * 80)
    print("Common mask summary")
    print(f"Total pixels: {common_mask.size}")
    print(f"Common valid pixels: {common_mask.sum()}")
    print(f"Fraction valid: {common_mask.sum() / common_mask.size:.4f}")

    # Save mask
    ref_header = headers[500].copy()
    ref_header["BUNIT"] = "1"
    ref_header["COMMENT"] = "Common valid mask: finite and positive in all four Herschel bands."

    mask_path = OUT_DIR / "common_valid_mask.fits"
    fits.writeto(mask_path, common_mask.astype(np.uint8), ref_header, overwrite=True)
    print(f"Saved mask FITS: {mask_path}")

    save_mask_plot(common_mask, ref_header)

    # Apply mask to all maps and save
    for wavelength, data in maps.items():
        clean = data.copy()
        clean[~common_mask] = np.nan

        header = headers[wavelength].copy()
        header["BUNIT"] = "MJy/sr"
        header["WAVELNTH"] = wavelength
        header["COMMENT"] = "Common mask applied: invalid or non-positive pixels set to NaN."

        out_fits = OUT_DIR / f"herschel_{wavelength}um_common.fits"
        fits.writeto(out_fits, clean, header, overwrite=True)

        print("=" * 80)
        print(f"Saved cleaned map: {out_fits}")
        print(f"Band: {wavelength} micron")
        print(f"Finite pixels after mask: {np.isfinite(clean).sum()} / {clean.size}")
        print(f"Min after mask: {np.nanmin(clean):.6e}")
        print(f"Max after mask: {np.nanmax(clean):.6e}")
        print(f"Median after mask: {np.nanmedian(clean):.6e}")

        save_map_plot(clean, header, wavelength)

    print("\nDone.")
    print(f"Common masked FITS files saved in: {OUT_DIR}")
    print(f"Plots saved in: {PLOT_DIR}")


if __name__ == "__main__":
    main()
