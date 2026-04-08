#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent.parent
PHOT_DIR = BASE_DIR / "output" / "photometry"
PLOT_DIR = BASE_DIR / "output" / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def plot_single_night(csv_path: Path, night_label: str):
    df = pd.read_csv(csv_path)

    if "JD" not in df.columns:
        raise ValueError(f"JD column not found in {csv_path}")

    jd0 = df["JD"].min()
    df["hours_from_start"] = (df["JD"] - jd0) * 24.0

    # Differential magnitude
    plt.figure(figsize=(8, 5))
    plt.plot(df["hours_from_start"], df["target_minus_ref_ensemble"], marker="o", linestyle="-")
    plt.xlabel("Hours from first exposure")
    plt.ylabel("Differential instrumental magnitude")
    plt.title(f"BL Lac Differential Light Curve ({night_label})")
    plt.tight_layout()
    out1 = PLOT_DIR / f"BLLac_diffmag_hours_{night_label}.png"
    plt.savefig(out1, dpi=200)
    plt.close()

    # Calibrated R magnitude
    plt.figure(figsize=(8, 5))
    plt.plot(df["hours_from_start"], df["target_Rmag"], marker="o", linestyle="-")
    plt.gca().invert_yaxis()
    plt.xlabel("Hours from first exposure")
    plt.ylabel("Target R magnitude")
    plt.title(f"BL Lac Light Curve ({night_label})")
    plt.tight_layout()
    out2 = PLOT_DIR / f"BLLac_Rmag_hours_{night_label}.png"
    plt.savefig(out2, dpi=200)
    plt.close()

    print(f"Saved:\n{out1}\n{out2}")


def plot_combined(csv_path: Path, label: str):
    df = pd.read_csv(csv_path)

    if "JD" not in df.columns:
        raise ValueError(f"JD column not found in {csv_path}")

    jd0 = df["JD"].min()
    df["hours_from_start"] = (df["JD"] - jd0) * 24.0

    plt.figure(figsize=(9, 5))
    plt.plot(df["hours_from_start"], df["target_minus_ref_ensemble"], marker="o", linestyle="-")
    plt.xlabel("Hours from first observation")
    plt.ylabel("Differential instrumental magnitude")
    plt.title(f"BL Lac Differential Light Curve ({label})")
    plt.tight_layout()
    out1 = PLOT_DIR / f"BLLac_diffmag_hours_{label}.png"
    plt.savefig(out1, dpi=200)
    plt.close()

    plt.figure(figsize=(9, 5))
    plt.plot(df["hours_from_start"], df["target_Rmag"], marker="o", linestyle="-")
    plt.gca().invert_yaxis()
    plt.xlabel("Hours from first observation")
    plt.ylabel("Target R magnitude")
    plt.title(f"BL Lac Light Curve ({label})")
    plt.tight_layout()
    out2 = PLOT_DIR / f"BLLac_Rmag_hours_{label}.png"
    plt.savefig(out2, dpi=200)
    plt.close()

    print(f"Saved:\n{out1}\n{out2}")


def main():
    parser = argparse.ArgumentParser(description="Make nicer light-curve plots.")
    parser.add_argument("--mode", choices=["single", "combined"], required=True)
    parser.add_argument("--label", required=True, help="Night label, e.g. 20200824 or 20200824-20200825")
    args = parser.parse_args()

    csv_path = PHOT_DIR / f"BLLac_{args.label}.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    if args.mode == "single":
        plot_single_night(csv_path, args.label)
    else:
        plot_combined(csv_path, args.label)


if __name__ == "__main__":
    main()