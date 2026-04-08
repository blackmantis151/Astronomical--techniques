from pathlib import Path
import numpy as np
from astropy.io import fits


# Speed of light in km/s
C_KMS = 299792.458


def load_sdss_spectrum(fits_path):
    """
    Load one SDSS spectrum from a FITS file.

    Parameters
    ----------
    fits_path : str or Path
        Path to FITS file

    Returns
    -------
    wavelength : np.ndarray
        Wavelength array in Angstrom
    flux : np.ndarray
        Flux array
    header : fits.Header
        Primary header
    """
    fits_path = Path(fits_path)

    with fits.open(fits_path) as hdul:
        header = hdul[0].header
        data = hdul[1].data

        # SDSS spectra usually store log10(lambda)
        loglam = data["loglam"]
        flux = data["flux"]

        wavelength = 10 ** loglam

    return wavelength, flux, header


def get_flux_stats(flux):
    """
    Basic statistics for a flux array.
    """
    flux = np.asarray(flux)
    return {
        "min": float(np.min(flux)),
        "max": float(np.max(flux)),
        "mean": float(np.mean(flux)),
        "median": float(np.median(flux)),
        "std": float(np.std(flux)),
    }


def find_nearest_index(array, value):
    """
    Return index of element in array closest to value.
    """
    array = np.asarray(array)
    return int(np.argmin(np.abs(array - value)))


def is_wavelength_in_range(wavelength, test_value):
    """
    Check whether a wavelength value lies inside the spectrum range.
    """
    wavelength = np.asarray(wavelength)
    return float(np.min(wavelength)) <= float(test_value) <= float(np.max(wavelength))


def moving_average(y, window=25):
    """
    Simple moving average smoothing.

    Parameters
    ----------
    y : array-like
        Input array
    window : int
        Smoothing window size

    Returns
    -------
    smooth_y : np.ndarray
        Smoothed array
    """
    y = np.asarray(y)

    if window < 1:
        raise ValueError("window must be >= 1")

    if window == 1:
        return y.copy()

    kernel = np.ones(window, dtype=float) / window
    return np.convolve(y, kernel, mode="same")


def compute_residual(flux, smooth_flux):
    """
    Residual = flux - smooth continuum estimate
    """
    flux = np.asarray(flux)
    smooth_flux = np.asarray(smooth_flux)
    return flux - smooth_flux


def cut_region(wavelength, flux, center, window=100):
    """
    Cut a local spectral region around a given wavelength center.

    Parameters
    ----------
    wavelength : array-like
    flux : array-like
    center : float
        Central wavelength in Angstrom
    window : float
        Half-width of region in Angstrom

    Returns
    -------
    w_cut : np.ndarray
    f_cut : np.ndarray
    """
    wavelength = np.asarray(wavelength)
    flux = np.asarray(flux)

    mask = (wavelength >= center - window) & (wavelength <= center + window)
    return wavelength[mask], flux[mask]


def estimate_local_continuum(w_cut, f_cut, edge_fraction=0.2):
    """
    Estimate local continuum using the edge pixels of a cut region.

    Idea:
    We assume the line sits near the middle, and the edges better represent
    the local baseline.

    Parameters
    ----------
    w_cut : np.ndarray
    f_cut : np.ndarray
    edge_fraction : float
        Fraction of points from each edge used to estimate continuum

    Returns
    -------
    continuum_level : float
    """
    w_cut = np.asarray(w_cut)
    f_cut = np.asarray(f_cut)

    n = len(f_cut)
    if n == 0:
        return None

    edge_n = max(1, int(edge_fraction * n))

    left = f_cut[:edge_n]
    right = f_cut[-edge_n:]

    continuum_level = float(np.median(np.concatenate([left, right])))
    return continuum_level


def subtract_local_continuum(w_cut, f_cut, edge_fraction=0.2):
    """
    Subtract local continuum from a cut spectral region.

    Returns
    -------
    f_sub : np.ndarray
        Continuum-subtracted flux
    continuum_level : float
        Estimated continuum level
    """
    continuum_level = estimate_local_continuum(w_cut, f_cut, edge_fraction=edge_fraction)
    if continuum_level is None:
        return None, None

    f_sub = np.asarray(f_cut) - continuum_level
    return f_sub, continuum_level


def estimate_line_center(wavelength, flux, expected_center, window=30, mode="emission"):
    """
    Rough first estimate of a line center from a local region.

    Parameters
    ----------
    wavelength : np.ndarray
    flux : np.ndarray
    expected_center : float
        Rough expected observed wavelength
    window : float
        Half-width of local search region
    mode : str
        "emission" -> choose local maximum
        "absorption" -> choose local minimum

    Returns
    -------
    obs_lambda : float or None
    obs_flux : float or None
    local_data : tuple or None
        (w_cut, f_cut)
    """
    w_cut, f_cut = cut_region(wavelength, flux, expected_center, window=window)

    if len(w_cut) == 0:
        return None, None, None

    if mode == "emission":
        idx_local = np.argmax(f_cut)
    elif mode == "absorption":
        idx_local = np.argmin(f_cut)
    else:
        raise ValueError("mode must be 'emission' or 'absorption'")

    obs_lambda = float(w_cut[idx_local])
    obs_flux = float(f_cut[idx_local])

    return obs_lambda, obs_flux, (w_cut, f_cut)


def refine_line_center_centroid(wavelength, flux, expected_center, window=20, mode="emission"):
    """
    Refine line center using a weighted centroid in a local region.

    This is still simple, but better than picking just one pixel.

    For emission:
        weights = positive continuum-subtracted flux
    For absorption:
        weights = positive depth below continuum

    Parameters
    ----------
    wavelength : np.ndarray
    flux : np.ndarray
    expected_center : float
    window : float
    mode : str
        "emission" or "absorption"

    Returns
    -------
    refined_lambda : float or None
    feature_strength : float or None
    local_data : tuple or None
        (w_cut, f_cut, f_feature, continuum_level)
    """
    w_cut, f_cut = cut_region(wavelength, flux, expected_center, window=window)

    if len(w_cut) < 3:
        return None, None, None

    f_sub, continuum_level = subtract_local_continuum(w_cut, f_cut)

    if f_sub is None:
        return None, None, None

    if mode == "emission":
        weights = np.clip(f_sub, 0, None)
    elif mode == "absorption":
        weights = np.clip(-f_sub, 0, None)
    else:
        raise ValueError("mode must be 'emission' or 'absorption'")

    total_weight = np.sum(weights)
    if total_weight <= 0:
        return None, None, None

    refined_lambda = float(np.sum(w_cut * weights) / total_weight)
    feature_strength = float(np.max(weights))

    return refined_lambda, feature_strength, (w_cut, f_cut, weights, continuum_level)


def compute_redshift(obs_lambda, rest_lambda):
    """
    Compute redshift from observed and rest wavelengths.

    z = (lambda_obs - lambda_rest) / lambda_rest
    """
    return (obs_lambda - rest_lambda) / rest_lambda


def compute_velocity_from_z(z):
    """
    Non-relativistic approximation: v = c z

    Good enough for this project unless your instructor explicitly asks for
    a relativistic correction.
    """
    return C_KMS * z


def compute_velocity_from_line(obs_lambda, rest_lambda):
    """
    Compute redshift and velocity from one spectral line.

    Returns
    -------
    z : float
    v : float
        Velocity in km/s
    """
    z = compute_redshift(obs_lambda, rest_lambda)
    v = compute_velocity_from_z(z)
    return z, v


def predict_observed_wavelength(rest_lambda, z):
    """
    Predict where another rest line should appear for the same redshift.

    lambda_obs = lambda_rest * (1 + z)
    """
    return rest_lambda * (1 + z)


def summarize_line_measurement(line_name, rest_lambda, obs_lambda):
    """
    Convenience helper for printing or logging line measurements.
    """
    z, v = compute_velocity_from_line(obs_lambda, rest_lambda)
    return {
        "line_name": line_name,
        "rest_lambda": float(rest_lambda),
        "obs_lambda": float(obs_lambda),
        "z": float(z),
        "velocity_kms": float(v),
    }


def weighted_linear_fit(x, y, yerr):
    """
    Weighted linear fit: y = m*x + b using numpy.polyfit.

    Parameters
    ----------
    x : array-like
    y : array-like
    yerr : array-like
        Uncertainties in y

    Returns
    -------
    m, b, m_err, b_err
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)

    if np.any(yerr <= 0):
        raise ValueError("All yerr values must be positive for weighted fit.")

    weights = 1.0 / yerr
    coeffs, cov = np.polyfit(x, y, deg=1, w=weights, cov=True)

    m, b = coeffs
    m_err = float(np.sqrt(cov[0, 0]))
    b_err = float(np.sqrt(cov[1, 1]))

    return float(m), float(b), m_err, b_err