from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import curve_fit


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_DIR = PROJECT_ROOT / "data" / "processed" / "common_masked"
OUT_DIR = PROJECT_ROOT / "data" / "processed" / "sed_fitting"
PLOT_DIR = PROJECT_ROOT / "plots" / "sed_fitting"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


MAP_FILES = {
    160: INPUT_DIR / "herschel_160um_common.fits",
    250: INPUT_DIR / "herschel_250um_common.fits",
    350: INPUT_DIR / "herschel_350um_common.fits",
    500: INPUT_DIR / "herschel_500um_common.fits",
}

MASK_FILE = INPUT_DIR / "common_valid_mask.fits"


# Physical constants in CGS
H = 6.62607015e-27       # erg s
K_B = 1.380649e-16       # erg K^-1
C = 2.99792458e10        # cm s^-1
M_H = 1.6735575e-24      # g

MU_H2 = 2.8              # mean molecular weight per H2 molecule
BETA = 2.0

# Dust opacity model
# kappa_nu = kappa0 * (nu / nu0)^beta
# kappa0 here is gas+dust opacity in cm^2/g
KAPPA0 = 0.1             # cm^2/g
NU0 = 1.0e12             # Hz, 1000 GHz


def planck_nu_cgs(nu_hz, temperature_k):
    """
    Planck function B_nu(T) in CGS units:
    erg s^-1 cm^-2 Hz^-1 sr^-1
    """
    x = H * nu_hz / (K_B * temperature_k)

    # Avoid overflow for very small/large trial values.
    x = np.clip(x, 1e-8, 700)

    return (2.0 * H * nu_hz**3 / C**2) / (np.exp(x) - 1.0)


def kappa_nu(nu_hz):
    """
    Dust opacity in cm^2/g.
    """
    return KAPPA0 * (nu_hz / NU0) ** BETA


def modified_blackbody_model(wavelength_um, temperature_k, log10_nh2):
    """
    Return intensity in MJy/sr.

    Input:
        wavelength_um: wavelength in micron
        temperature_k: dust temperature in K
        log10_nh2: log10 of H2 column density in cm^-2

    Model:
        I_nu = B_nu(T) * kappa_nu * mu_H2 * m_H * N(H2)

    Units:
        B_nu is in erg s^-1 cm^-2 Hz^-1 sr^-1.
        1 MJy/sr = 1e-17 erg s^-1 cm^-2 Hz^-1 sr^-1.
    """
    wavelength_cm = wavelength_um * 1e-4
    nu_hz = C / wavelength_cm

    nh2 = 10.0 ** log10_nh2
    sigma = MU_H2 * M_H * nh2

    intensity_cgs = planck_nu_cgs(nu_hz, temperature_k) * kappa_nu(nu_hz) * sigma

    intensity_mjysr = intensity_cgs / 1e-17

    return intensity_mjysr


def read_map(path):
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(float)
        header = hdul[0].header.copy()
    return data, header


def robust_limits(image):
    finite = np.isfinite(image)

    if finite.sum() == 0:
        return 0, 1

    values = image[finite]
    vmin = np.nanpercentile(values, 2)
    vmax = np.nanpercentile(values, 98)

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin, vmax = np.nanmin(values), np.nanmax(values)

    return vmin, vmax


def save_fits(data, header, filename, bunit, comment):
    out_header = header.copy()
    out_header["BUNIT"] = bunit
    out_header["COMMENT"] = comment
    fits.writeto(OUT_DIR / filename, data, out_header, overwrite=True)


def save_map_plot(data, header, title, cbar_label, filename, cmap="inferno"):
    fig = plt.figure(figsize=(7, 6))

    try:
        wcs = WCS(header)
        ax = plt.subplot(projection=wcs)
        im = ax.imshow(data, origin="lower", cmap=cmap)
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    except Exception:
        ax = plt.subplot()
        im = ax.imshow(data, origin="lower", cmap=cmap)
        ax.set_xlabel("Pixel x")
        ax.set_ylabel("Pixel y")

    vmin, vmax = robust_limits(data)
    im.set_clim(vmin, vmax)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(PLOT_DIR / filename, dpi=220)
    plt.close()

    print(f"Saved plot: {PLOT_DIR / filename}")


def save_example_sed_plot(wavelengths_um, observed, fitted, temperature, nh2, y, x):
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.scatter(wavelengths_um, observed, label="Observed", s=60)
    ax.plot(wavelengths_um, fitted, label="Best fit", linewidth=2)

    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Surface brightness [MJy/sr]")
    ax.set_title(
        f"Example SED fit at pixel (y={y}, x={x})\n"
        f"T = {temperature:.2f} K, N(H2) = {nh2:.2e} cm$^{{-2}}$"
    )

    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    out = PLOT_DIR / "example_sed_fit.png"
    plt.savefig(out, dpi=220)
    plt.close()

    print(f"Saved example SED plot: {out}")


def main():
    print("\nFitting modified blackbody SED pixel by pixel\n")

    wavelengths_um = np.array([160.0, 250.0, 350.0, 500.0])

    maps = {}
    headers = {}

    for wave in [160, 250, 350, 500]:
        path = MAP_FILES[wave]

        if not path.exists():
            raise FileNotFoundError(f"Missing input map: {path}")

        data, header = read_map(path)
        maps[wave] = data
        headers[wave] = header

        print("=" * 80)
        print(f"Read {wave} micron map")
        print(f"Shape: {data.shape}")
        print(f"Finite pixels: {np.isfinite(data).sum()} / {data.size}")
        print(f"Median: {np.nanmedian(data):.6e} MJy/sr")

    reference_header = headers[500].copy()
    shape = maps[500].shape

    cube = np.stack([maps[160], maps[250], maps[350], maps[500]], axis=0)

    valid = np.all(np.isfinite(cube), axis=0)
    valid &= np.all(cube > 0, axis=0)

    print("=" * 80)
    print("SED fitting mask")
    print(f"Total pixels: {valid.size}")
    print(f"Valid pixels for fitting: {valid.sum()}")

    temperature_map = np.full(shape, np.nan, dtype=float)
    nh2_map = np.full(shape, np.nan, dtype=float)
    log_nh2_map = np.full(shape, np.nan, dtype=float)
    chi2_map = np.full(shape, np.nan, dtype=float)

    # Bounds:
    # temperature from 5 K to 60 K
    # log10 N(H2) from 18 to 25
    lower_bounds = [5.0, 18.0]
    upper_bounds = [60.0, 25.0]

    p0 = [15.0, 21.5]

    n_success = 0
    n_fail = 0

    ys, xs = np.where(valid)

    for count, (y, x) in enumerate(zip(ys, xs), start=1):
        obs = cube[:, y, x]

        try:
            popt, pcov = curve_fit(
                modified_blackbody_model,
                wavelengths_um,
                obs,
                p0=p0,
                bounds=(lower_bounds, upper_bounds),
                maxfev=10000,
            )

            temp, log_nh2 = popt
            model = modified_blackbody_model(wavelengths_um, temp, log_nh2)

            residual = obs - model
            chi2 = np.nanmean((residual / obs) ** 2)

            temperature_map[y, x] = temp
            log_nh2_map[y, x] = log_nh2
            nh2_map[y, x] = 10.0 ** log_nh2
            chi2_map[y, x] = chi2

            n_success += 1

        except Exception:
            n_fail += 1

        if count % 1000 == 0:
            print(f"Processed {count} / {len(ys)} valid pixels")

    print("=" * 80)
    print("Fit summary")
    print(f"Successful fits: {n_success}")
    print(f"Failed fits: {n_fail}")

    print("\nTemperature map:")
    print(f"Min T: {np.nanmin(temperature_map):.3f} K")
    print(f"Max T: {np.nanmax(temperature_map):.3f} K")
    print(f"Median T: {np.nanmedian(temperature_map):.3f} K")

    print("\nColumn density map:")
    print(f"Min N(H2): {np.nanmin(nh2_map):.3e} cm^-2")
    print(f"Max N(H2): {np.nanmax(nh2_map):.3e} cm^-2")
    print(f"Median N(H2): {np.nanmedian(nh2_map):.3e} cm^-2")

    save_fits(
        temperature_map,
        reference_header,
        "dust_temperature_map.fits",
        "K",
        "Dust temperature from modified blackbody fit. beta fixed to 2."
    )

    save_fits(
        nh2_map,
        reference_header,
        "column_density_NH2_map.fits",
        "cm^-2",
        "Molecular hydrogen column density from modified blackbody fit. beta fixed to 2."
    )

    save_fits(
        log_nh2_map,
        reference_header,
        "log_column_density_NH2_map.fits",
        "log10(cm^-2)",
        "log10 molecular hydrogen column density from modified blackbody fit."
    )

    save_fits(
        chi2_map,
        reference_header,
        "relative_chi2_map.fits",
        "1",
        "Mean squared relative residual of modified blackbody fit."
    )

    save_map_plot(
        temperature_map,
        reference_header,
        "Dust temperature map",
        "Dust temperature [K]",
        "dust_temperature_map.png",
        cmap="inferno",
    )

    save_map_plot(
        np.log10(nh2_map),
        reference_header,
        "Column density map",
        r"log$_{10}$ N(H$_2$) [cm$^{-2}$]",
        "column_density_logNH2_map.png",
        cmap="viridis",
    )

    save_map_plot(
        chi2_map,
        reference_header,
        "Relative fit residual map",
        "Mean squared relative residual",
        "relative_chi2_map.png",
        cmap="magma",
    )

    # Make one example SED plot from the brightest valid 500 micron pixel.
    candidate = np.where(valid, maps[500], np.nan)

    if np.isfinite(candidate).sum() > 0:
        y0, x0 = np.unravel_index(np.nanargmax(candidate), candidate.shape)

        observed = cube[:, y0, x0]
        temp = temperature_map[y0, x0]
        nh2 = nh2_map[y0, x0]
        fitted = modified_blackbody_model(wavelengths_um, temp, np.log10(nh2))

        save_example_sed_plot(
            wavelengths_um,
            observed,
            fitted,
            temp,
            nh2,
            y0,
            x0,
        )

    print("\nDone.")
    print(f"SED fitting outputs saved in: {OUT_DIR}")
    print(f"Plots saved in: {PLOT_DIR}")


if __name__ == "__main__":
    main()
