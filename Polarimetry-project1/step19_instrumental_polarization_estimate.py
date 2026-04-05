import numpy as np
import pandas as pd

from config import POL_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories


def estimate_instrumental_polarization():
    ensure_directories(OUTPUT_SUBDIRS)

    df = pd.read_csv(POL_DIR / "polarimetry_results.csv")

    q_mean = df["q"].mean()
    u_mean = df["u"].mean()

    p_inst = np.sqrt(q_mean**2 + u_mean**2)
    theta_inst = 0.5 * np.degrees(np.arctan2(u_mean, q_mean))

    out_txt = POL_DIR / "instrumental_polarization_estimate.txt"
    with open(out_txt, "w") as f:
        f.write("Crude field-average instrumental polarization estimate\n")
        f.write("=====================================================\n")
        f.write(f"Mean q = {q_mean}\n")
        f.write(f"Mean u = {u_mean}\n")
        f.write(f"P_inst = {p_inst}\n")
        f.write(f"theta_inst_deg = {theta_inst}\n")

    print(f"Saved: {out_txt}")
    print(f"Mean q = {q_mean:.6f}")
    print(f"Mean u = {u_mean:.6f}")
    print(f"P_inst = {p_inst:.6f}")
    print(f"theta_inst = {theta_inst:.3f} deg")


if __name__ == "__main__":
    estimate_instrumental_polarization()