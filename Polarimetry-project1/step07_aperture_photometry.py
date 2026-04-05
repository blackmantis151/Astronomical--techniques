import pandas as pd
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry

from config import (
    STACK_DIR, POL_DIR, PHOT_DIR,
    APERTURE_RADII, ANNULUS_R_IN, ANNULUS_R_OUT,
    OUTPUT_SUBDIRS
)
from utils import ensure_directories, load_fits_data


def measure_flux(data, x, y, r_ap, r_in, r_out):
    ap = CircularAperture((x, y), r=r_ap)
    ann = CircularAnnulus((x, y), r_in=r_in, r_out=r_out)

    ap_sum = aperture_photometry(data, ap)["aperture_sum"][0]

    ann_mask = ann.to_mask(method="center")
    ann_data = ann_mask.multiply(data)
    ann_vals = ann_data[ann_mask.data > 0]

    sky = float(pd.Series(ann_vals).median())
    flux = ap_sum - sky * ap.area
    return flux, sky


def run_photometry():
    ensure_directories(OUTPUT_SUBDIRS)

    #pairs = pd.read_csv(POL_DIR / "beam_pairs.csv")
    pairs = pd.read_csv(POL_DIR / "beam_pairs_filtered.csv")  # Use filtered pairs for photometry

    rows = []

    for g in [1, 2, 3, 4]:
        data, _ = load_fits_data(STACK_DIR / f"group_{g}_stack.fits")

        for _, p in pairs.iterrows():
            for r_ap in APERTURE_RADII:
                f_o, sky_o = measure_flux(data, p["x_o"], p["y_o"], r_ap, ANNULUS_R_IN, ANNULUS_R_OUT)
                f_e, sky_e = measure_flux(data, p["x_e"], p["y_e"], r_ap, ANNULUS_R_IN, ANNULUS_R_OUT)

                rows.append({
                    "pair_id": int(p["pair_id"]),
                    "group_id": g,
                    "r_ap": r_ap,
                    "x_o": p["x_o"],
                    "y_o": p["y_o"],
                    "x_e": p["x_e"],
                    "y_e": p["y_e"],
                    "flux_o": f_o,
                    "flux_e": f_e,
                    "sky_o": sky_o,
                    "sky_e": sky_e,
                })

    df = pd.DataFrame(rows)
    out_csv = PHOT_DIR / "aperture_photometry.csv"
    df.to_csv(out_csv, index=False)

    print(f"Saved photometry table: {out_csv}")


if __name__ == "__main__":
    run_photometry()