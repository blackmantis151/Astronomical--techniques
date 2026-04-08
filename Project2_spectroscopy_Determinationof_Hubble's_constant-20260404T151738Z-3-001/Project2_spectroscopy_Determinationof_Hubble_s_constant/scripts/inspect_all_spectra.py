from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from functions import load_sdss_spectrum
from config import GALAXIES

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_hubble"
OUT_DIR = BASE_DIR / "output"
PLOT_DIR = OUT_DIR / "plots"
TEXT_DIR = OUT_DIR / "text"

PLOT_DIR.mkdir(parents=True, exist_ok=True)
TEXT_DIR.mkdir(parents=True, exist_ok=True)


def moving_average(y, window=25):
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


def analyze_spectrum(wavelength, flux):
    stats = {
        "min": float(np.min(flux)),
        "max": float(np.max(flux)),
        "mean": float(np.mean(flux)),
        "median": float(np.median(flux)),
        "std": float(np.std(flux)),
    }

    smooth_flux = moving_average(flux, window=25)
    residual = flux - smooth_flux

    peak_residual = float(np.max(residual))
    dip_residual = float(np.min(residual))

    return stats, smooth_flux, residual, peak_residual, dip_residual


def classify_spectrum(stats, peak_residual, dip_residual):
    std = stats["std"]

    flags = []
    if peak_residual > 3 * std:
        flags.append("strong_emission_features")
    if abs(dip_residual) > 3 * std:
        flags.append("strong_absorption_features")
    if peak_residual < 2 * std and abs(dip_residual) < 2 * std:
        flags.append("mostly_smooth_or_weak_lines")

    if not flags:
        flags.append("mixed_or_unclear")

    return flags


def plot_inspection(filename, wavelength, flux, smooth_flux, residual, ref_rest_lambda, ref_name):
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)

    axes[0].plot(wavelength, flux, lw=0.7)
    axes[0].axvline(ref_rest_lambda, ls="--", lw=1, label=f"{ref_name} rest λ = {ref_rest_lambda:.1f} Å")
    axes[0].set_title(f"{filename} : raw spectrum")
    axes[0].set_ylabel("Flux")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(wavelength, smooth_flux, lw=1.0)
    axes[1].axvline(ref_rest_lambda, ls="--", lw=1)
    axes[1].set_title("Smooth continuum estimate")
    axes[1].set_ylabel("Flux")
    axes[1].grid(alpha=0.3)

    axes[2].plot(wavelength, residual, lw=0.7)
    axes[2].axhline(0, color="black", lw=0.8)
    axes[2].axvline(ref_rest_lambda, ls="--", lw=1)
    axes[2].set_title("Residual = flux - smooth continuum")
    axes[2].set_xlabel("Wavelength (Angstrom)")
    axes[2].set_ylabel("Residual")
    axes[2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOT_DIR / f"{filename}_inspection.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    lines = []
    lines.append("INSPECTION REPORT FOR ALL SPECTRA\n")
    lines.append("=" * 80 + "\n")

    for gal in GALAXIES:
        filename = gal["filename"]
        ref_name = gal["reference_line"]
        ref_rest_lambda = gal["rest_lambda"]
        distance_mpc = gal["distance_mpc"]
        distance_err = gal["distance_err"]

        filepath = DATA_DIR / filename
        wavelength, flux, header = load_sdss_spectrum(filepath)

        stats, smooth_flux, residual, peak_residual, dip_residual = analyze_spectrum(wavelength, flux)
        flags = classify_spectrum(stats, peak_residual, dip_residual)

        plot_inspection(filename, wavelength, flux, smooth_flux, residual, ref_rest_lambda, ref_name)

        lines.append(f"FILE: {filename}\n")
        lines.append(f"Distance: {distance_mpc:.2f} +/- {distance_err:.2f} Mpc\n")
        lines.append(f"Reference line: {ref_name}\n")
        lines.append(f"Reference rest wavelength: {ref_rest_lambda:.1f} Å\n")
        lines.append(f"Wavelength range: {wavelength.min():.2f} to {wavelength.max():.2f} Å\n")
        lines.append(f"Number of points: {len(wavelength)}\n")
        lines.append(f"Flux min: {stats['min']:.4f}\n")
        lines.append(f"Flux max: {stats['max']:.4f}\n")
        lines.append(f"Flux mean: {stats['mean']:.4f}\n")
        lines.append(f"Flux median: {stats['median']:.4f}\n")
        lines.append(f"Flux std: {stats['std']:.4f}\n")
        lines.append(f"Max residual: {peak_residual:.4f}\n")
        lines.append(f"Min residual: {dip_residual:.4f}\n")
        lines.append(f"Classification: {', '.join(flags)}\n")
        lines.append("Manual task: inspect the plot and write a rough observed wavelength for the shifted reference line in config.py\n")
        lines.append(f"Saved plot: output/plots/{filename}_inspection.png\n")
        lines.append("-" * 80 + "\n")

    report_path = TEXT_DIR / "inspect_all_spectra_report.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    print(f"Saved report to: {report_path}")


if __name__ == "__main__":
    main()