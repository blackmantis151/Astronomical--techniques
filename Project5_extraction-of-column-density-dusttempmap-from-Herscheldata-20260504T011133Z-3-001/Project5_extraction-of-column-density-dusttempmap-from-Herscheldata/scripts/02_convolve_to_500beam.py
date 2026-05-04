from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.convolution import convolve_fft


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_DIR = PROJECT_ROOT / "data" / "processed" / "extracted_mjysr"

KERNEL_DIR = (
    PROJECT_ROOT
    / "data"
    / "Dust_temp_demo-20260504T115037Z-3-001"
    / "Dust_temp_demo"
    / "Kernels"
)

OUT_DIR = PROJECT_ROOT / "data" / "processed" / "convolved_to_500"
PLOT_DIR = PROJECT_ROOT / "plots" / "convolved_maps"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


MAP_FILES = {
    160: INPUT_DIR / "herschel_160um_mjysr.fits",
    250: INPUT_DIR / "herschel_250um_mjysr.fits",
    350: INPUT_DIR / "herschel_350um_mjysr.fits",
    500: INPUT_DIR / "herschel_500um_mjysr.fits",
}


KERNEL_FILES = {
    160: KERNEL_DIR / "Kernel_LoRes_PACS_160_to_SPIRE_500_regrided_to_PACS_160.fits",
    250: KERNEL_DIR / "Kernel_LoRes_SPIRE_250_to_SPIRE_500.fits",
    350: KERNEL_DIR / "Kernel_LoRes_SPIRE_350_to_SPIRE_500.fits",
}


def read_fits_image(path):
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(float)
        header = hdul[0].header.copy()
    return data, header


def read_kernel(path):
    with fits.open(path) as hdul:
        kernel = hdul[0].data.astype(float)

    # Normalize kernel so that total flux/surface brightness scale is preserved.
    kernel_sum = np.nansum(kernel)

    if not np.isfinite(kernel_sum) or abs(kernel_sum) < 1e-12:
        raise ValueError(f"Kernel sum is bad for {path.name}: {kernel_sum}")

    kernel = kernel / kernel_sum

    return kernel


def convolve_image(data, kernel):
    """
    Convolve image with provided beam-matching kernel.

    preserve_nan=True keeps invalid regions as NaN as much as possible.
    nan_treatment='interpolate' avoids NaNs spreading through the full image.
    normalize_kernel=False because we already normalized the kernel manually.
    """
    convolved = convolve_fft(
        data,
        kernel,
        boundary="fill",
        fill_value=np.nan,
        nan_treatment="interpolate",
        normalize_kernel=False,
        preserve_nan=True,
        allow_huge=True,
    )

    return convolved


def robust_limits(image):
    finite = np.isfinite(image)
    positive = finite & (image > 0)

    if positive.sum() > 100:
        values = image[positive]
    else:
        values = image[finite]

    vmin = np.nanpercentile(values, 2)
    vmax = np.nanpercentile(values, 99)

    return vmin, vmax


def save_plot(image, wavelength):
    fig, ax = plt.subplots(figsize=(7, 6))

    im = ax.imshow(image, origin="lower", cmap="inferno")
    vmin, vmax = robust_limits(image)
    im.set_clim(vmin, vmax)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Surface brightness [MJy/sr]")

    ax.set_xlabel("Pixel x")
    ax.set_ylabel("Pixel y")
    ax.set_title(f"{wavelength} micron map convolved to 500 micron beam")

    out_png = PLOT_DIR / f"convolved_{wavelength}um_to_500beam.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"Saved plot: {out_png}")


def main():
    print("\nConvolving Herschel maps to common 500 micron beam\n")

    for wavelength, map_path in MAP_FILES.items():
        print("=" * 80)
        print(f"Band: {wavelength} micron")
        print(f"Input map: {map_path}")

        if not map_path.exists():
            raise FileNotFoundError(f"Missing input map: {map_path}")

        data, header = read_fits_image(map_path)

        print(f"Input shape: {data.shape}")
        print(f"Input unit: {header.get('BUNIT')}")
        print(f"Input finite pixels: {np.isfinite(data).sum()} / {data.size}")
        print(f"Input min: {np.nanmin(data):.6e}")
        print(f"Input max: {np.nanmax(data):.6e}")
        print(f"Input median: {np.nanmedian(data):.6e}")

        if wavelength == 500:
            print("500 micron map is already the target beam. Copying unchanged.")
            output_data = data.copy()
            header["COMMENT"] = "500 micron map used as common target beam."
        else:
            kernel_path = KERNEL_FILES[wavelength]

            if not kernel_path.exists():
                raise FileNotFoundError(f"Missing kernel: {kernel_path}")

            kernel = read_kernel(kernel_path)

            print(f"Kernel: {kernel_path.name}")
            print(f"Kernel shape: {kernel.shape}")
            print(f"Kernel sum after normalization: {np.nansum(kernel):.6e}")
            print(f"Kernel min: {np.nanmin(kernel):.6e}")
            print(f"Kernel max: {np.nanmax(kernel):.6e}")

            output_data = convolve_image(data, kernel)

            header["COMMENT"] = f"Convolved from {wavelength} micron native beam to SPIRE 500 micron beam."
            header["COMMENT"] = f"Kernel used: {kernel_path.name}"

        header["BUNIT"] = "MJy/sr"
        header["WAVELNTH"] = wavelength
        header["BEAMMAT"] = "SPIRE500"

        out_fits = OUT_DIR / f"herschel_{wavelength}um_mjysr_conv500.fits"
        fits.writeto(out_fits, output_data, header, overwrite=True)

        print(f"Output file: {out_fits}")
        print(f"Output finite pixels: {np.isfinite(output_data).sum()} / {output_data.size}")
        print(f"Output min: {np.nanmin(output_data):.6e}")
        print(f"Output max: {np.nanmax(output_data):.6e}")
        print(f"Output median: {np.nanmedian(output_data):.6e}")

        save_plot(output_data, wavelength)

    print("\nDone.")
    print(f"Convolved FITS files saved in: {OUT_DIR}")
    print(f"Plots saved in: {PLOT_DIR}")


if __name__ == "__main__":
    main()
