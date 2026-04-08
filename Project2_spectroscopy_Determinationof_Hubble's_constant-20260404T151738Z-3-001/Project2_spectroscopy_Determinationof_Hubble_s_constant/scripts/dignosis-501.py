from pathlib import Path
import matplotlib.pyplot as plt

from functions import load_sdss_spectrum
from config import REST_LINES

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_hubble"

filename = "spec-0501-52235-0386.fits"
filepath = DATA_DIR / filename

# Accept manual distance and cleaned H0
distance_mpc = 1597.44
H0 = 67.3299
c_kms = 299792.458

v = H0 * distance_mpc
z = v / c_kms

predicted = {
    "Hbeta": REST_LINES["Hbeta"] * (1 + z),
    "OIII_5007": REST_LINES["OIII_5007"] * (1 + z),
    "Halpha": REST_LINES["Halpha"] * (1 + z),
    "SII_6733": REST_LINES["SII_6733"] * (1 + z),
}

wavelength, flux, header = load_sdss_spectrum(filepath)

plt.figure(figsize=(13, 5))
plt.plot(wavelength, flux, lw=0.7)

for name, lam in predicted.items():
    plt.axvline(lam, ls="--", lw=1.2, label=f"{name} pred = {lam:.1f} Å")

plt.xlim(6200, 9300)
plt.xlabel("Wavelength (Angstrom)")
plt.ylabel("Flux")
plt.title(f"{filename} diagnostic using manual distance + cleaned H0")
plt.grid(alpha=0.3)
plt.legend()
plt.show()

print("v =", v, "km/s")
print("z =", z)
for k, val in predicted.items():
    print(k, val)