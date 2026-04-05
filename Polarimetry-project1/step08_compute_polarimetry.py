import numpy as np
import pandas as pd

from config import GROUP_TO_HWP, PHOT_DIR, POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def compute_qu(f0_o, f0_e, f45_o, f45_e, f22_o, f22_e, f67_o, f67_e):
    vals = [f0_o, f0_e, f45_o, f45_e, f22_o, f22_e, f67_o, f67_e]

    if any((not np.isfinite(v)) or (v <= 0) for v in vals):
        return np.nan, np.nan, "non_positive_or_invalid_flux"

    rq_arg = (f0_o / f0_e) / (f45_o / f45_e)
    ru_arg = (f22_o / f22_e) / (f67_o / f67_e)

    if not np.isfinite(rq_arg) or not np.isfinite(ru_arg):
        return np.nan, np.nan, "invalid_ratio_argument"

    if rq_arg <= 0:
        return np.nan, np.nan, "non_positive_rq_argument"

    if ru_arg <= 0:
        return np.nan, np.nan, "non_positive_ru_argument"

    rq = np.sqrt(rq_arg)
    q = (rq - 1) / (rq + 1)

    ru = np.sqrt(ru_arg)
    u = (ru - 1) / (ru + 1)

    if (not np.isfinite(q)) or (not np.isfinite(u)):
        return np.nan, np.nan, "invalid_q_or_u"

    return q, u, None


def compute_p_theta(q, u):
    if not np.isfinite(q) or not np.isfinite(u):
        return np.nan, np.nan, "invalid_q_or_u"

    p = np.sqrt(q**2 + u**2)
    theta = 0.5 * np.degrees(np.arctan2(u, q))

    if not np.isfinite(p) or not np.isfinite(theta):
        return np.nan, np.nan, "invalid_p_or_theta"

    return p, theta, None


def run_polarimetry(chosen_aperture=5):
    ensure_directories(OUTPUT_SUBDIRS)

    df = pd.read_csv(PHOT_DIR / "aperture_photometry.csv")
    df = df[df["r_ap"] == chosen_aperture].copy()

    hwp_map = {g: GROUP_TO_HWP[g] for g in GROUP_TO_HWP}
    df["hwp"] = df["group_id"].map(hwp_map)

    required = [0.0, 22.5, 45.0, 67.5]

    valid_results = []
    rejected_results = []

    for pid in sorted(df["pair_id"].unique()):
        s = df[df["pair_id"] == pid].copy()

        angle_map = {}
        for _, row in s.iterrows():
            angle_map[row["hwp"]] = row

        if not all(a in angle_map for a in required):
            missing = [a for a in required if a not in angle_map]
            rejected_results.append({
                "pair_id": pid,
                "reason": f"missing_hwp_angles_{missing}",
                "r_ap": chosen_aperture
            })
            continue

        f0_o  = angle_map[0.0]["flux_o"]
        f0_e  = angle_map[0.0]["flux_e"]
        f45_o = angle_map[45.0]["flux_o"]
        f45_e = angle_map[45.0]["flux_e"]
        f22_o = angle_map[22.5]["flux_o"]
        f22_e = angle_map[22.5]["flux_e"]
        f67_o = angle_map[67.5]["flux_o"]
        f67_e = angle_map[67.5]["flux_e"]

        q, u, qu_reason = compute_qu(
            f0_o, f0_e,
            f45_o, f45_e,
            f22_o, f22_e,
            f67_o, f67_e,
        )

        if qu_reason is not None:
            rejected_results.append({
                "pair_id": pid,
                "reason": qu_reason,
                "r_ap": chosen_aperture,
                "f0_o": f0_o,
                "f0_e": f0_e,
                "f45_o": f45_o,
                "f45_e": f45_e,
                "f22_o": f22_o,
                "f22_e": f22_e,
                "f67_o": f67_o,
                "f67_e": f67_e,
            })
            continue

        p, theta, pt_reason = compute_p_theta(q, u)

        if pt_reason is not None:
            rejected_results.append({
                "pair_id": pid,
                "reason": pt_reason,
                "r_ap": chosen_aperture,
                "q": q,
                "u": u
            })
            continue

        valid_results.append({
            "pair_id": pid,
            "x": angle_map[0.0]["x_o"],
            "y": angle_map[0.0]["y_o"],
            "q": q,
            "u": u,
            "P": p,
            "theta_deg": theta,
            "r_ap": chosen_aperture,
        })

    valid_df = pd.DataFrame(valid_results)
    rejected_df = pd.DataFrame(rejected_results)

    # final safety cleaning for valid results
    if len(valid_df) > 0:
        valid_df = valid_df.replace([np.inf, -np.inf], np.nan)
        valid_df = valid_df.dropna(subset=["q", "u", "P", "theta_deg"]).copy()

    valid_csv = POL_DIR / "polarimetry_results.csv"
    rejected_csv = POL_DIR / "polarimetry_rejected_pairs.csv"

    valid_df.to_csv(valid_csv, index=False)
    rejected_df.to_csv(rejected_csv, index=False)

    print(f"Saved valid polarimetry results: {valid_csv}")
    print(f"Saved rejected pair log: {rejected_csv}")
    print(f"Valid stars: {len(valid_df)}")
    print(f"Rejected stars: {len(rejected_df)}")

    if len(valid_df) > 0:
        print("\nFirst few valid results:")
        print(valid_df.head())

    if len(rejected_df) > 0:
        print("\nFirst few rejected results:")
        print(rejected_df.head())


if __name__ == "__main__":
    run_polarimetry()