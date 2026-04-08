from pathlib import Path
import matplotlib.pyplot as plt

from functions import load_sdss_spectrum
from config import GALAXIES

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_hubble"

INDEX = 2  # change this manually

gal = GALAXIES[INDEX]
filename = gal["filename"]
ref_name = gal["reference_line"]
rest_lambda = gal["rest_lambda"]
guess_obs_lambda = gal["guess_obs_lambda"]

filepath = DATA_DIR / filename
wavelength, flux, header = load_sdss_spectrum(filepath)

plt.figure(figsize=(13, 5))
plt.plot(wavelength, flux, lw=0.7)
plt.axvline(rest_lambda, ls="--", lw=1.2, label=f"{ref_name} rest λ = {rest_lambda:.1f} Å")

if guess_obs_lambda is not None:
    plt.axvline(guess_obs_lambda, ls=":", lw=1.2, label=f"manual guess = {guess_obs_lambda:.1f} Å")

plt.xlabel("Wavelength (Angstrom)")
plt.ylabel("Flux")
plt.title(f"{filename} | reference line: {ref_name}")
plt.grid(alpha=0.3)
plt.legend()
plt.show()

print(f"File: {filename}")
print(f"Reference line: {ref_name}")
print(f"Rest wavelength: {rest_lambda:.1f} Å")
print("Inspect the full plot.")
print("Choose a plausible shifted feature and write it into config.py as guess_obs_lambda.")