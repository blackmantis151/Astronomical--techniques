from pathlib import Path
import matplotlib.pyplot as plt


def plot_full_spectrum(wavelength, flux, title="", savepath=None, show=True):
    """
    Plot the full spectrum.
    """
    plt.figure(figsize=(12, 5))
    plt.plot(wavelength, flux, lw=0.8)
    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel("Flux")
    plt.title(title)
    plt.grid(alpha=0.3)

    if savepath is not None:
        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(savepath, dpi=200, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()


def plot_line_region(w_cut, f_cut, expected_center=None, detected_center=None,
                     title="", savepath=None, show=True):
    """
    Plot a zoomed in spectral region around one line.
    """
    plt.figure(figsize=(9, 4))
    plt.plot(w_cut, f_cut, lw=1.0)

    if expected_center is not None:
        plt.axvline(expected_center, ls="--", lw=1, label=f"Expected: {expected_center:.2f}")

    if detected_center is not None:
        plt.axvline(detected_center, ls=":", lw=1.2, label=f"Detected: {detected_center:.2f}")

    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel("Flux")
    plt.title(title)
    plt.grid(alpha=0.3)
    plt.legend()

    if savepath is not None:
        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(savepath, dpi=200, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()


def plot_hubble_diagram(distance, velocity, velocity_err=None,
                        slope=None, intercept=None,
                        title="Hubble Diagram", savepath=None, show=True):
    """
    Plot velocity vs distance with optional best fit line.
    """
    plt.figure(figsize=(8, 6))

    if velocity_err is not None:
        plt.errorbar(distance, velocity, yerr=velocity_err, fmt="o", capsize=4, label="Galaxies")
    else:
        plt.plot(distance, velocity, "o", label="Galaxies")

    if slope is not None and intercept is not None:
        xfit = [min(distance), max(distance)]
        yfit = [slope * x + intercept for x in xfit]
        plt.plot(xfit, yfit, label=f"Fit: v = {slope:.2f} d + {intercept:.2f}")

    plt.xlabel("Distance (Mpc)")
    plt.ylabel("Velocity (km/s)")
    plt.title(title)
    plt.grid(alpha=0.3)
    plt.legend()

    if savepath is not None:
        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(savepath, dpi=200, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()