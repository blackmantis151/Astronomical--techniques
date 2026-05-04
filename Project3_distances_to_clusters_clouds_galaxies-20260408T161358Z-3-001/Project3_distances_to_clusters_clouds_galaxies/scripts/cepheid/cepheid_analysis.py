from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

CEP_DIR = ROOT / "cepheid"

catalog_file = CEP_DIR / "cepF.dat"
i_file = CEP_DIR / "OGLE-LMC-CEP-0017.dat"
v_file = CEP_DIR / "OGLE-LMC-CEP-0017.dat.1"

cepheid_id = "OGLE-LMC-CEP-0017"


def read_lightcurve(path):
    return pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["time", "mag", "mag_err"]
    ).dropna()


def intensity_mean_mag(mags):
    flux = 10 ** (-0.4 * np.asarray(mags))
    return -2.5 * np.log10(np.mean(flux))


def read_ogle_catalog_line(catalog_path, target_id):
    """
    For OGLE LMC fundamental Cepheids cepF.dat, the line format is:

    ID  I  V  P  dP  T0  A_I  A_V  R21  phi21  ...
    
    For your object:
    OGLE-LMC-CEP-0017  15.237 15.899   3.6772904 ...
    
    Therefore:
    parts[1] = mean I
    parts[2] = mean V
    parts[3] = period
    """
    with open(catalog_path, "r") as f:
        for line in f:
            if line.strip().startswith(target_id):
                parts = line.split()
                return {
                    "id": parts[0],
                    "I_cat": float(parts[1]),
                    "V_cat": float(parts[2]),
                    "P": float(parts[3]),
                    "raw_line": line.strip()
                }

    raise ValueError(f"{target_id} not found in {catalog_path}")


# -----------------------------
# Read data
# -----------------------------
lc_I = read_lightcurve(i_file)
lc_V = read_lightcurve(v_file)

cat = read_ogle_catalog_line(catalog_file, cepheid_id)

P = cat["P"]
m_I_catalog = cat["I_cat"]
m_V_catalog = cat["V_cat"]

print("\nCatalog line found:")
print(cat["raw_line"])

print("\nParsed catalog values:")
print(f"ID = {cat['id']}")
print(f"Mean I = {m_I_catalog:.4f}")
print(f"Mean V = {m_V_catalog:.4f}")
print(f"Period P = {P:.7f} days")

print("\nI-band light curve:")
print(lc_I.head())

print("\nV-band light curve:")
print(lc_V.head())

# -----------------------------
# Mean magnitudes
# -----------------------------
mean_I_mag = lc_I["mag"].mean()
mean_V_mag = lc_V["mag"].mean()

mean_I_int = intensity_mean_mag(lc_I["mag"])
mean_V_int = intensity_mean_mag(lc_V["mag"])

print("\nMean magnitudes from light curves:")
print(f"I arithmetic mean = {mean_I_mag:.4f}")
print(f"I intensity mean  = {mean_I_int:.4f}")
print(f"V arithmetic mean = {mean_V_mag:.4f}")
print(f"V intensity mean  = {mean_V_int:.4f}")

# Use OGLE catalog mean V because catalog means are already standardized.
m_V = m_V_catalog

print(f"\nUsing catalog mean V magnitude m_V = {m_V:.4f}")

# -----------------------------
# Cepheid PL relation
# -----------------------------
M_V = -2.76 * np.log10(P) - 1.40

A_V = 0.30

mu = m_V - M_V - A_V
d_pc = 10 ** ((mu + 5) / 5)
d_kpc = d_pc / 1000

mu_lmc_ref = 18.49
d_lmc_ref_kpc = 49.9

print("\nCepheid distance result:")
print(f"M_V = {M_V:.4f}")
print(f"A_V = {A_V:.2f} mag")
print(f"Distance modulus mu = {mu:.4f} mag")
print(f"Distance = {d_kpc:.3f} kpc")
print(f"Reference LMC distance = {d_lmc_ref_kpc:.1f} kpc")
print(f"Difference = {d_kpc - d_lmc_ref_kpc:.2f} kpc")

# -----------------------------
# Plot 1: time light curve
# -----------------------------
plt.figure(figsize=(8, 5))
plt.errorbar(lc_I["time"], lc_I["mag"], yerr=lc_I["mag_err"],
             fmt=".", alpha=0.7, label="I band")
plt.errorbar(lc_V["time"], lc_V["mag"], yerr=lc_V["mag_err"],
             fmt=".", alpha=0.9, label="V band")
plt.gca().invert_yaxis()
plt.xlabel("Time")
plt.ylabel("Magnitude")
plt.title(f"{cepheid_id}: OGLE Light Curve")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "cepheid_time_lightcurve.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 2: phase folded light curve
# -----------------------------
lc_I["phase"] = (lc_I["time"] / P) % 1
lc_V["phase"] = (lc_V["time"] / P) % 1

plt.figure(figsize=(8, 5))
plt.errorbar(lc_I["phase"], lc_I["mag"], yerr=lc_I["mag_err"],
             fmt=".", alpha=0.7, label="I band")
plt.errorbar(lc_V["phase"], lc_V["mag"], yerr=lc_V["mag_err"],
             fmt=".", alpha=0.9, label="V band")

plt.errorbar(lc_I["phase"] + 1, lc_I["mag"], yerr=lc_I["mag_err"],
             fmt=".", alpha=0.35)
plt.errorbar(lc_V["phase"] + 1, lc_V["mag"], yerr=lc_V["mag_err"],
             fmt=".", alpha=0.45)

plt.gca().invert_yaxis()
plt.xlabel("Phase")
plt.ylabel("Magnitude")
plt.title(f"{cepheid_id}: Phase-folded Light Curve, P = {P:.3f} d")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "cepheid_phase_folded.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 3: PL relation point
# -----------------------------
logP = np.log10(P)
logP_grid = np.linspace(0.0, 2.0, 200)
M_grid = -2.76 * logP_grid - 1.40

plt.figure(figsize=(7, 6))
plt.plot(logP_grid, M_grid, label=r"$M_V=-2.76\log P - 1.40$")
plt.scatter([logP], [M_V], s=80, label=cepheid_id)
plt.gca().invert_yaxis()
plt.xlabel(r"$\log_{10}(P/\mathrm{day})$")
plt.ylabel(r"$M_V$")
plt.title("Cepheid Period Luminosity Relation")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "cepheid_PL_relation_point.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 4: extinction sensitivity
# -----------------------------
A_values = np.linspace(0.0, 0.8, 100)
mu_values = m_V - M_V - A_values
d_values_kpc = 10 ** ((mu_values + 5) / 5) / 1000

plt.figure(figsize=(7, 5))
plt.plot(A_values, d_values_kpc)
plt.axvline(A_V, linestyle="--", label=f"Adopted A_V = {A_V:.2f}")
plt.axhline(d_lmc_ref_kpc, linestyle=":", label="Reference LMC distance")
plt.xlabel(r"$A_V$ [mag]")
plt.ylabel("Distance [kpc]")
plt.title("Effect of V-band Extinction on Cepheid Distance")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "cepheid_extinction_sensitivity.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# Save result
# -----------------------------
result_file = CEP_DIR / "cepheid_distance_result.txt"

with open(result_file, "w") as f:
    f.write(f"Cepheid ID: {cepheid_id}\n")
    f.write(f"Catalog line: {cat['raw_line']}\n")
    f.write(f"Period P = {P:.7f} days\n")
    f.write(f"Catalog mean I = {m_I_catalog:.4f}\n")
    f.write(f"Catalog mean V used = {m_V:.4f}\n")
    f.write(f"PL relation: M_V = -2.76 log10(P) - 1.40\n")
    f.write(f"M_V = {M_V:.4f}\n")
    f.write(f"A_V = {A_V:.2f} mag\n")
    f.write(f"Distance modulus mu = {mu:.4f} mag\n")
    f.write(f"Distance = {d_kpc:.3f} kpc\n")
    f.write(f"Reference LMC distance = {d_lmc_ref_kpc:.1f} kpc\n")
    f.write(f"Difference = {d_kpc - d_lmc_ref_kpc:.2f} kpc\n")

print(f"\nSaved result file: {result_file}")
print(f"Saved plots in: {PLOT_DIR}")