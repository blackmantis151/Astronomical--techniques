from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS


PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = (
    PROJECT_ROOT
    / "data"
    / "Dust_temp_demo-20260504T115037Z-3-001"
    / "Dust_temp_demo"
    / "Original_files"
)

OUT_DIR = PROJECT_ROOT / "data" / "processed" / "extracted_mjysr"
PLOT_DIR = PROJECT_ROOT / "plots" / "raw_maps"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


# Approximate Herschel beam FWHM values in arcsec.
# These are used only for converting Jy/beam to MJy/sr.
# PACS 160 is Jy/pixel, so it uses pixel solid angle instead.
BEAM_FWHM_ARCSEC = {
    250: 18.2,
    350: 24.9,
    500: 36.3,
}


def wavelength_from_filename_or_header(filename, primary_header):
    name = filename.lower()

    if "pmp250" in name:
        return 250
    if "pmp350" in name:
        return 350
    if "pmp500" in name:
        return 500
    if "pacs" in name or "hppjsmapr" in name:
        return 160

    wave = primary_header.get("WAVELNTH")
    if wave is None:
        raise ValueError(f"Could not identify wavelength for {filename}")

    return int(round(float(wave)))


def pixel_area_sr(header):
    """
    Compute pixel solid angle from FITS WCS header.

    CDELT1 and CDELT2 are in degrees per pixel.
    Pixel area in steradian = |CDELT1*CDELT2| * (pi/180)^2
    """
    cdelt1 = header.get("CDELT1")
    cdelt2 = header.get("CDELT2")

    if cdelt1 is None or cdelt2 is None:
        raise ValueError("CDELT1/CDELT2 missing. Cannot compute pixel area.")

    return abs(float(cdelt1) * float(cdelt2)) * (np.pi / 180.0) ** 2


def beam_area_sr_from_fwhm(fwhm_arcsec):
    """
    Gaussian beam solid angle.

    Omega_beam = 1.1331 * FWHM^2

    FWHM must be in radians.
    """
    fwhm_rad = fwhm_arcsec / 206265.0
    return 1.1331 * fwhm_rad ** 2


def convert_to_mjysr(data, header, bunit, wavelength):
    """
    Convert the image to MJy/sr.

    Jy/pixel:
        I[Jy/sr] = image[Jy/pixel] / pixel_area[sr]

    Jy/beam:
        I[Jy/sr] = image[Jy/beam] / beam_area[sr]

    Then:
        MJy/sr = Jy/sr / 1e6
    """
    if bunit is None:
        raise ValueError("BUNIT missing.")

    bunit_clean = bunit.strip().lower().replace(" ", "")

    if bunit_clean in ["mjy/sr", "mjy/sr-1", "mjy/sr"]:
        return data.copy()

    if bunit_clean in ["jy/pixel", "jy/pix"]:
        omega_pix = pixel_area_sr(header)
        return data / omega_pix / 1e6

    if bunit_clean in ["jy/beam"]:
        if wavelength not in BEAM_FWHM_ARCSEC:
            raise ValueError(f"No beam FWHM defined for wavelength {wavelength} micron.")

        omega_beam = beam_area_sr_from_fwhm(BEAM_FWHM_ARCSEC[wavelength])
        return data / omega_beam / 1e6

    raise ValueError(f"Unknown unit: {bunit}")


def robust_vmin_vmax(image):
    finite = np.isfinite(image)
    positive = finite & (image > 0)

    if positive.sum() > 100:
        values = image[positive]
    else:
        values = image[finite]

    vmin = np.nanpercentile(values, 2)
    vmax = np.nanpercentile(values, 99)

    return vmin, vmax


def save_preview(image_mjysr, header, wavelength):
    fig = plt.figure(figsize=(7, 6))

    try:
        wcs = WCS(header)
        ax = plt.subplot(projection=wcs)
        im = ax.imshow(image_mjysr, origin="lower", cmap="inferno")
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    except Exception:
        ax = plt.subplot()
        im = ax.imshow(image_mjysr, origin="lower", cmap="inferno")
        ax.set_xlabel("Pixel x")
        ax.set_ylabel("Pixel y")

    vmin, vmax = robust_vmin_vmax(image_mjysr)
    im.set_clim(vmin, vmax)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Surface brightness [MJy/sr]")

    ax.set_title(f"Herschel {wavelength} micron map")

    out_png = PLOT_DIR / f"raw_{wavelength}um_mjysr.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"Saved preview plot: {out_png}")


def main():
    files = sorted(RAW_DIR.glob("*.fits"))

    if len(files) == 0:
        raise FileNotFoundError(f"No FITS files found in {RAW_DIR}")

    print("\nExtracting science images and converting to MJy/sr\n")

    for file in files:
        print("=" * 80)
        print(f"Input file: {file.name}")

        with fits.open(file) as hdul:
            primary_header = hdul[0].header

            image_hdu = hdul["image"]
            data = image_hdu.data.astype(float)
            header = image_hdu.header.copy()

            wavelength = wavelength_from_filename_or_header(file.name, primary_header)
            bunit = header.get("BUNIT")

            print(f"Wavelength: {wavelength} micron")
            print(f"Original unit: {bunit}")
            print(f"Original shape: {data.shape}")

            image_mjysr = convert_to_mjysr(data, header, bunit, wavelength)

            header["BUNIT"] = "MJy/sr"
            header["WAVELNTH"] = wavelength
            header["COMMENT"] = "Science image extracted from HDU image"
            header["COMMENT"] = "Converted to MJy/sr for SED fitting"

            out_fits = OUT_DIR / f"herschel_{wavelength}um_mjysr.fits"

            fits.writeto(out_fits, image_mjysr, header, overwrite=True)

            finite = np.isfinite(image_mjysr)
            print(f"Output file: {out_fits}")
            print(f"Finite pixels: {finite.sum()} / {image_mjysr.size}")
            print(f"Min MJy/sr: {np.nanmin(image_mjysr):.6e}")
            print(f"Max MJy/sr: {np.nanmax(image_mjysr):.6e}")
            print(f"Median MJy/sr: {np.nanmedian(image_mjysr):.6e}")

            save_preview(image_mjysr, header, wavelength)

    print("\nDone.")
    print(f"Converted files saved in: {OUT_DIR}")
    print(f"Preview plots saved in: {PLOT_DIR}")


if __name__ == "__main__":
    main()
