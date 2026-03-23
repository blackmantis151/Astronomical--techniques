import numpy as np
import pandas as pd

from config import GROUP_TO_HWP, PHOT_DIR, POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def compute_qu(f0_o, f0_e, f45_o, f45_e, f22_o, f22_e, f67_o, f67_e):
    rq = np.sqrt((f0_o / f0_e) / (f45_o / f45_e))
    q = (rq - 1) / (rq + 1)

    ru = np.sqrt((f22_o / f22_e) / (f67_o / f67_e))
    u = (ru - 1) / (ru + 1)

    return q, u


def compute_p_theta(q, u):
    p = np.sqrt(q**2 + u**2)
    theta = 0.5 * np.degrees(np.arctan2(u, q))
    return p, theta


def run_polarimetry(chosen_aperture=5):
    ensure_directories(OUTPUT_SUBDIRS)

    df = pd.read_csv(PHOT_DIR / "aperture_photometry.csv")
    df = df[df["r_ap"] == chosen_aperture].copy()

    hwp_map = {g: GROUP_TO_HWP[g] for g in GROUP_TO_HWP}
    df["hwp"] = df["group_id"].map(hwp_map)

    results = []

    for pid in sorted(df["pair_id"].unique()):
        s = df[df["pair_id"] == pid]

        angle_map = {}
        for _, row in s.iterrows():
            angle_map[row["hwp"]] = row

        required = [0.0, 22.5, 45.0, 67.5]
        if not all(a in angle_map for a in required):
            continue

        q, u = compute_qu(
            angle_map[0.0]["flux_o"],   angle_map[0.0]["flux_e"],
            angle_map[45.0]["flux_o"],  angle_map[45.0]["flux_e"],
            angle_map[22.5]["flux_o"],  angle_map[22.5]["flux_e"],
            angle_map[67.5]["flux_o"],  angle_map[67.5]["flux_e"],
        )
        p, theta = compute_p_theta(q, u)

        results.append({
            "pair_id": pid,
            "x": angle_map[0.0]["x_o"],
            "y": angle_map[0.0]["y_o"],
            "q": q,
            "u": u,
            "P": p,
            "theta_deg": theta,
            "r_ap": chosen_aperture,
        })

    out_df = pd.DataFrame(results)
    out_csv = POL_DIR / "polarimetry_results.csv"
    out_df.to_csv(out_csv, index=False)

    print(f"Saved polarimetry results: {out_csv}")
    print(out_df.head())


if __name__ == "__main__":
    run_polarimetry()