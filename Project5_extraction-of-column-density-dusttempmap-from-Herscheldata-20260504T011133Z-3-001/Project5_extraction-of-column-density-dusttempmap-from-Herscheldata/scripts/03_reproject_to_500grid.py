from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_DIR = PROJECT_ROOT / "data" / "processed" / "convolved_to_500"
OUT_DIR = PROJECT_ROOT / "data" / "processed" / "reprojected_to_500grid"
PLOT_DIR = PROJECT_ROOT / "plots" / "reprojected_maps"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


MAP_FILES = {
    160: INPUT_DIR / "herschel_160um_mjysr_conv500.fits",
    250: INPUT_DIR / "herschel_250um_mjysr_conv500.fits",
    350: INPUT_DIR / "herschel_350um_mjysr_conv500.fits",
    500: INPUT_DIR / "herschel_500um_mjysr_conv500.fits",
}


REFERENCE_WAVELENGTH = 500
REFERENCE_FILE = MAP_FILES[REFERENCE_WAVELENGTH]


def read_fits(path):
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
        vmin, vmax = np.nanmin(values), np.nanmax(values)

    return vmin, vmax


def save_plot(image, header, wavelength):
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

    ax.set_title(f"{wavelength} micron reprojected to 500 micron grid")

    out_png = PLOT_DIR / f"reprojected_{wavelength}um_to_500grid.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"Saved plot: {out_png}")


def main():
    print("\nReprojecting all convolved maps to the 500 micron grid\n")

    if not REFERENCE_FILE.exists():
        raise FileNotFoundError(f"Reference file missing: {REFERENCE_FILE}")

    ref_data, ref_header = read_fits(REFERENCE_FILE)
    ref_wcs = WCS(ref_header)
    ref_shape = ref_data.shape

    print("=" * 80)
    print("Reference grid")
    print(f"Reference file: {REFERENCE_FILE}")
    print(f"Reference shape: {ref_shape}")
    print(f"Reference CDELT1: {ref_header.get('CDELT1')}")
    print(f"Reference CDELT2: {ref_header.get('CDELT2')}")
    print(f"Reference BUNIT: {ref_header.get('BUNIT')}")

    for wavelength, path in MAP_FILES.items():
        print("=" * 80)
        print(f"Band: {wavelength} micron")
        print(f"Input file: {path}")

        if not path.exists():
            raise FileNotFoundError(f"Missing input file: {path}")

        data, header = read_fits(path)

        print(f"Input shape: {data.shape}")
        print(f"Input finite pixels: {np.isfinite(data).sum()} / {data.size}")

        if wavelength == REFERENCE_WAVELENGTH:
            print("500 micron is the reference grid. Copying unchanged.")
            reproj_data = data.copy()
            reproj_header = ref_header.copy()
            footprint = np.isfinite(reproj_data).astype(float)
        else:
            reproj_data, footprint = reproject_interp(
                (data, WCS(header)),
                ref_header,
                shape_out=ref_shape,
                order="bilinear",
            )

            reproj_header = ref_header.copy()
            reproj_header["BUNIT"] = "MJy/sr"
            reproj_header["WAVELNTH"] = wavelength
            reproj_header["COMMENT"] = f"Reprojected from {wavelength} micron native grid to 500 micron grid."

        # Remove pixels that were not covered by reprojection.
        reproj_data[footprint <= 0] = np.nan

        out_fits = OUT_DIR / f"herschel_{wavelength}um_mjysr_conv500_regrid500.fits"
        fits.writeto(out_fits, reproj_data, reproj_header, overwrite=True)

        print(f"Output file: {out_fits}")
        print(f"Output shape: {reproj_data.shape}")
        print(f"Output finite pixels: {np.isfinite(reproj_data).sum()} / {reproj_data.size}")
        print(f"Output min: {np.nanmin(reproj_data):.6e}")
        print(f"Output max: {np.nanmax(reproj_data):.6e}")
        print(f"Output median: {np.nanmedian(reproj_data):.6e}")

        save_plot(reproj_data, reproj_header, wavelength)

    print("\nDone.")
    print(f"Reprojected FITS files saved in: {OUT_DIR}")
    print(f"Plots saved in: {PLOT_DIR}")


if __name__ == "__main__":
    main()
