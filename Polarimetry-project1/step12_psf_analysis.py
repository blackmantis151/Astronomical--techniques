import numpy as np
import matplotlib.pyplot as plt

from astropy.visualization import simple_norm

from config import (
    STACK_DIR,
    PLOT_DIR,
    OUTPUT_SUBDIRS,
    ALIGNMENT_GROUP_ID,
    ALIGNMENT_X,
    ALIGNMENT_Y,
    ALIGNMENT_CUTOUT_HALF_SIZE,
)
from utils import ensure_directories, load_fits_data, extract_cutout


def radial_profile(data, x0, y0):
    """
    Compute azimuthally averaged radial profile around (x0, y0).
    """
    y, x = np.indices(data.shape)
    r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    r_int = r.astype(int)

    tbin = np.bincount(r_int.ravel(), data.ravel())
    nr = np.bincount(r_int.ravel())

    profile = tbin / nr
    return profile


def compute_fwhm(profile):
    """
    Approximate FWHM from radial profile.
    Returns diameter in pixels.
    """
    peak = np.max(profile)
    half_max = peak / 2.0

    below = np.where(profile < half_max)[0]
    if len(below) == 0:
        return None

    r_half = below[0]
    return 2.0 * r_half


def second_moment_ellipticity(data):
    """
    Estimate ellipticity from second intensity moments.
    Returns ellipticity, major_sigma, minor_sigma.
    """
    img = np.array(data, dtype=float)
    img = img - np.median(img)
    img[img < 0] = 0

    total = img.sum()
    if total <= 0:
        return np.nan, np.nan, np.nan

    y, x = np.indices(img.shape)
    x0 = (x * img).sum() / total
    y0 = (y * img).sum() / total

    x2 = (((x - x0) ** 2) * img).sum() / total
    y2 = (((y - y0) ** 2) * img).sum() / total
    xy = (((x - x0) * (y - y0)) * img).sum() / total

    cov = np.array([[x2, xy], [xy, y2]])
    eigvals = np.linalg.eigvalsh(cov)

    major = np.sqrt(np.max(eigvals))
    minor = np.sqrt(np.min(eigvals))

    if major <= 0:
        return np.nan, np.nan, np.nan

    ellipticity = 1.0 - (minor / major)
    return ellipticity, major, minor


def run_psf_analysis():
    ensure_directories(OUTPUT_SUBDIRS)

    stack_path = STACK_DIR / f"group_{ALIGNMENT_GROUP_ID}_stack.fits"
    data, _ = load_fits_data(stack_path)

    cutout, x1, y1 = extract_cutout(
        data,
        ALIGNMENT_X,
        ALIGNMENT_Y,
        half_size=ALIGNMENT_CUTOUT_HALF_SIZE,
    )

    x_local = ALIGNMENT_X - x1
    y_local = ALIGNMENT_Y - y1

    # 1. Save cutout image
    norm = simple_norm(cutout, stretch="linear", percent=99.5)
    plt.figure(figsize=(6, 6))
    plt.imshow(cutout, origin="lower", cmap="gray", norm=norm)
    plt.plot(x_local, y_local, marker="+", color="yellow", markersize=12, mew=2)
    plt.xlabel("X pixel in cutout")
    plt.ylabel("Y pixel in cutout")
    plt.title("PSF cutout of alignment star")
    out_cut = PLOT_DIR / "psf_cutout.png"
    plt.savefig(out_cut, dpi=150, bbox_inches="tight")
    plt.close()

    # 2. Save contour plot
    plt.figure(figsize=(6, 6))
    plt.imshow(cutout, origin="lower", cmap="gray", norm=norm, alpha=0.8)
    plt.contour(cutout, colors="cyan", linewidths=1)
    plt.plot(x_local, y_local, marker="+", color="yellow", markersize=12, mew=2)
    plt.xlabel("X pixel in cutout")
    plt.ylabel("Y pixel in cutout")
    plt.title("PSF contour plot")
    out_contour = PLOT_DIR / "psf_contour.png"
    plt.savefig(out_contour, dpi=150, bbox_inches="tight")
    plt.close()

    # 3. Radial profile
    profile = radial_profile(cutout, x_local, y_local)
    fwhm = compute_fwhm(profile)

    plt.figure(figsize=(7, 5))
    plt.plot(profile, marker="o")
    plt.axhline(np.max(profile) / 2.0, color="red", linestyle="--", label="Half maximum")
    if fwhm is not None:
        plt.axvline(fwhm / 2.0, color="green", linestyle="--", label=f"FWHM/2 ≈ {fwhm/2:.1f} pix")
    plt.xlabel("Radius (pixels)")
    plt.ylabel("Mean intensity")
    plt.title(f"Radial profile (FWHM ≈ {fwhm:.2f} pix)" if fwhm is not None else "Radial profile")
    plt.legend()
    plt.tight_layout()
    out_profile = PLOT_DIR / "psf_radial_profile.png"
    plt.savefig(out_profile, dpi=150, bbox_inches="tight")
    plt.close()

    # 4. Ellipticity
    ellipticity, major, minor = second_moment_ellipticity(cutout)

    # 5. Save summary text
    out_txt = PLOT_DIR / "psf_summary.txt"
    with open(out_txt, "w") as f:
        f.write("PSF analysis summary\n")
        f.write("====================\n")
        f.write(f"Group used: {ALIGNMENT_GROUP_ID}\n")
        f.write(f"Alignment star guess: x={ALIGNMENT_X}, y={ALIGNMENT_Y}\n")
        f.write(f"Cutout half-size: {ALIGNMENT_CUTOUT_HALF_SIZE}\n")
        f.write(f"Estimated FWHM (pixels): {fwhm}\n")
        f.write(f"Ellipticity: {ellipticity}\n")
        f.write(f"Major sigma-like width: {major}\n")
        f.write(f"Minor sigma-like width: {minor}\n")

    print(f"Saved: {out_cut}")
    print(f"Saved: {out_contour}")
    print(f"Saved: {out_profile}")
    print(f"Saved: {out_txt}")
    print(f"Estimated FWHM: {fwhm}")
    print(f"Estimated ellipticity: {ellipticity}")


if __name__ == "__main__":
    run_psf_analysis()