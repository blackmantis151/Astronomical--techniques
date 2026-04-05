import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.visualization import simple_norm
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.centroids import centroid_2dg

from config import STACK_DIR, PLOT_DIR, POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, load_fits_data


# =========================
# CONFIG: SET THESE
# =========================

# Initial guess for central O and E beam coordinates
# Update these if needed from your beam-offset plots
X_O_INIT, Y_O_INIT = 205, 216
DX_INIT, DY_INIT = 31, 5

X_E_INIT = X_O_INIT + DX_INIT
Y_E_INIT = Y_O_INIT + DY_INIT

# HWP mapping
HWP_MAP = {
    1: 67.5,
    2: 45.0,
    3: 22.5,
    4: 0.0
}

# Apertures to test
APERTURES = [5, 6, 8, 10, 12]

# Sky annulus
R_IN = 14
R_OUT = 22

# Centroid refinement cutout half-size
CENTROID_BOX = 12


# =========================
# HELPER FUNCTIONS
# =========================

def refine_centroid(data, x_init, y_init, box=CENTROID_BOX):
    """
    Refine centroid around an initial guess using 2D Gaussian centroiding.
    """
    x0 = int(round(x_init))
    y0 = int(round(y_init))

    x1 = max(0, x0 - box)
    x2 = min(data.shape[1], x0 + box + 1)
    y1 = max(0, y0 - box)
    y2 = min(data.shape[0], y0 + box + 1)

    cutout = data[y1:y2, x1:x2]

    yc, xc = centroid_2dg(cutout)

    x_refined = x1 + xc
    y_refined = y1 + yc

    return x_refined, y_refined, cutout


def do_aperture_photometry(data, x, y, r_ap, r_in=R_IN, r_out=R_OUT):
    """
    Circular aperture photometry with local sky subtraction from annulus.
    """
    aper = CircularAperture((x, y), r=r_ap)
    ann = CircularAnnulus((x, y), r_in=r_in, r_out=r_out)

    aper_sum = aperture_photometry(data, aper)["aperture_sum"][0]
    ann_sum = aperture_photometry(data, ann)["aperture_sum"][0]

    ann_area = ann.area
    aper_area = aper.area

    sky_mean = ann_sum / ann_area
    flux = aper_sum - sky_mean * aper_area

    return flux, sky_mean


def safe_compute_qu(f0_o, f0_e, f45_o, f45_e, f22_o, f22_e, f67_o, f67_e):
    vals = [f0_o, f0_e, f45_o, f45_e, f22_o, f22_e, f67_o, f67_e]

    if any((not np.isfinite(v)) or (v <= 0) for v in vals):
        return np.nan, np.nan

    rq_arg = (f0_o / f0_e) / (f45_o / f45_e)
    ru_arg = (f22_o / f22_e) / (f67_o / f67_e)

    if (not np.isfinite(rq_arg)) or (not np.isfinite(ru_arg)):
        return np.nan, np.nan

    if rq_arg <= 0 or ru_arg <= 0:
        return np.nan, np.nan

    rq = np.sqrt(rq_arg)
    ru = np.sqrt(ru_arg)

    q = (rq - 1) / (rq + 1)
    u = (ru - 1) / (ru + 1)

    return q, u


def safe_compute_p_theta(q, u):
    if not np.isfinite(q) or not np.isfinite(u):
        return np.nan, np.nan

    p = np.sqrt(q**2 + u**2)
    theta = 0.5 * np.degrees(np.arctan2(u, q))

    return p, theta


def load_all_groups():
    data_dict = {}
    for g in [1, 2, 3, 4]:
        path = STACK_DIR / f"group_{g}_stack.fits"
        data, _ = load_fits_data(path)
        data_dict[g] = data
    return data_dict


# =========================
# MAIN ANALYSIS
# =========================

def run_target_polarimetry():
    ensure_directories(OUTPUT_SUBDIRS)

    data_dict = load_all_groups()
    rows = []

    for r_ap in APERTURES:
        flux_dict = {}

        for g in [1, 2, 3, 4]:
            data = data_dict[g]

            # Refine centroids on this stacked image
            x_o_ref, y_o_ref, cutout_o = refine_centroid(data, X_O_INIT, Y_O_INIT)
            x_e_ref, y_e_ref, cutout_e = refine_centroid(data, X_E_INIT, Y_E_INIT)

            # Photometry using refined centroids
            f_o, sky_o = do_aperture_photometry(data, x_o_ref, y_o_ref, r_ap)
            f_e, sky_e = do_aperture_photometry(data, x_e_ref, y_e_ref, r_ap)

            angle = HWP_MAP[g]

            flux_dict[angle] = {
                "f_o": f_o,
                "f_e": f_e
            }

            rows.append({
                "group": g,
                "hwp": angle,
                "r_ap": r_ap,
                "x_o_ref": x_o_ref,
                "y_o_ref": y_o_ref,
                "x_e_ref": x_e_ref,
                "y_e_ref": y_e_ref,
                "flux_o": f_o,
                "flux_e": f_e,
                "sky_o": sky_o,
                "sky_e": sky_e
            })

        # Compute Stokes parameters for this aperture
        f0 = flux_dict[0.0]
        f45 = flux_dict[45.0]
        f22 = flux_dict[22.5]
        f67 = flux_dict[67.5]

        q, u = safe_compute_qu(
            f0["f_o"], f0["f_e"],
            f45["f_o"], f45["f_e"],
            f22["f_o"], f22["f_e"],
            f67["f_o"], f67["f_e"]
        )

        P, theta = safe_compute_p_theta(q, u)

        rows.append({
            "group": "RESULT",
            "hwp": -1,
            "r_ap": r_ap,
            "q": q,
            "u": u,
            "P": P,
            "theta": theta
        })

    df = pd.DataFrame(rows)

    out_csv = POL_DIR / "target_star_polarimetry.csv"
    df.to_csv(out_csv, index=False)

    print(f"Saved: {out_csv}")

    plot_diagnostics(df, data_dict)


# =========================
# PLOTS
# =========================

def plot_diagnostics(df, data_dict):
    result_df = df[df["group"] == "RESULT"].copy()

    # ---- P vs aperture ----
    plt.figure(figsize=(6, 4))
    plt.plot(result_df["r_ap"], result_df["P"], "o-")
    plt.xlabel("Aperture radius")
    plt.ylabel("P")
    plt.title("Target polarization vs aperture")
    plt.grid(True)
    out = PLOT_DIR / "target_P_vs_aperture.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")

    # ---- q vs aperture ----
    plt.figure(figsize=(6, 4))
    plt.plot(result_df["r_ap"], result_df["q"], "o-", label="q")
    plt.plot(result_df["r_ap"], result_df["u"], "s-", label="u")
    plt.xlabel("Aperture radius")
    plt.ylabel("Stokes parameter")
    plt.title("Target q and u vs aperture")
    plt.legend()
    plt.grid(True)
    out = PLOT_DIR / "target_qu_vs_aperture.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")

    # ---- Modulation curves ----
    for r_ap in APERTURES:
        sub = df[(df["r_ap"] == r_ap) & (df["group"] != "RESULT")].copy()
        sub = sub.sort_values("hwp")

        angles = sub["hwp"].values
        f_o = sub["flux_o"].values
        f_e = sub["flux_e"].values

        ratio = f_o / f_e
        diffnorm = (f_o - f_e) / (f_o + f_e)

        plt.figure(figsize=(6, 4))
        plt.plot(angles, f_o, "o-", label="O beam")
        plt.plot(angles, f_e, "s-", label="E beam")
        plt.xlabel("HWP angle (deg)")
        plt.ylabel("Flux")
        plt.title(f"Target flux vs HWP angle (r={r_ap})")
        plt.legend()
        plt.grid(True)
        out = PLOT_DIR / f"target_flux_vs_hwp_r{r_ap}.png"
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {out}")

        plt.figure(figsize=(6, 4))
        plt.plot(angles, ratio, "o-")
        plt.xlabel("HWP angle (deg)")
        plt.ylabel("O / E")
        plt.title(f"Target O/E ratio vs HWP angle (r={r_ap})")
        plt.grid(True)
        out = PLOT_DIR / f"target_ratio_vs_hwp_r{r_ap}.png"
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {out}")

        plt.figure(figsize=(6, 4))
        plt.plot(angles, diffnorm, "o-")
        plt.xlabel("HWP angle (deg)")
        plt.ylabel("(O - E) / (O + E)")
        plt.title(f"Target modulation vs HWP angle (r={r_ap})")
        plt.grid(True)
        out = PLOT_DIR / f"target_modulation_r{r_ap}.png"
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {out}")

    # ---- Visualization with refined centroids on group 1 ----
    data = data_dict[1]
    x_o_ref, y_o_ref, _ = refine_centroid(data, X_O_INIT, Y_O_INIT)
    x_e_ref, y_e_ref, _ = refine_centroid(data, X_E_INIT, Y_E_INIT)

    norm = simple_norm(data, stretch="linear", percent=99.5)
    plt.figure(figsize=(7, 7))
    plt.imshow(data, origin="lower", cmap="gray", norm=norm)
    plt.scatter([x_o_ref], [y_o_ref], color="yellow", label="O refined")
    plt.scatter([x_e_ref], [y_e_ref], color="cyan", label="E refined")
    plt.legend()
    plt.title("Target star beams (refined centroids)")
    out = PLOT_DIR / "target_beams.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# =========================
# RUN
# =========================

if __name__ == "__main__":
    run_target_polarimetry()