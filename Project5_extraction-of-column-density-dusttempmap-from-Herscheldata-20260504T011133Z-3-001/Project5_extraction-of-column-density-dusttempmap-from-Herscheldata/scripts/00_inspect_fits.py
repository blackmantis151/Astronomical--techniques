from pathlib import Path
from astropy.io import fits
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = (
    PROJECT_ROOT
    / "data"
    / "Dust_temp_demo-20260504T115037Z-3-001"
    / "Dust_temp_demo"
    / "Original_files"
)

KERNEL_DIR = (
    PROJECT_ROOT
    / "data"
    / "Dust_temp_demo-20260504T115037Z-3-001"
    / "Dust_temp_demo"
    / "Kernels"
)

fits_files = sorted(RAW_DIR.glob("*.fits"))
kernel_files = sorted(KERNEL_DIR.glob("*.fits"))


def inspect_hdu(hdu, index):
    header = hdu.header
    data = hdu.data

    print(f"\nHDU {index}: {hdu.name}")
    print(f"  Type    = {type(hdu).__name__}")
    print(f"  EXTNAME = {header.get('EXTNAME')}")
    print(f"  BUNIT   = {header.get('BUNIT')}")
    print(f"  NAXIS   = {header.get('NAXIS')}")
    print(f"  NAXIS1  = {header.get('NAXIS1')}")
    print(f"  NAXIS2  = {header.get('NAXIS2')}")
    print(f"  CTYPE1  = {header.get('CTYPE1')}")
    print(f"  CTYPE2  = {header.get('CTYPE2')}")
    print(f"  CDELT1  = {header.get('CDELT1')}")
    print(f"  CDELT2  = {header.get('CDELT2')}")
    print(f"  CRVAL1  = {header.get('CRVAL1')}")
    print(f"  CRVAL2  = {header.get('CRVAL2')}")
    print(f"  CRPIX1  = {header.get('CRPIX1')}")
    print(f"  CRPIX2  = {header.get('CRPIX2')}")
    print(f"  WAVELNTH= {header.get('WAVELNTH')}")
    print(f"  FILTER  = {header.get('FILTER')}")
    print(f"  BAND    = {header.get('BAND')}")
    print(f"  BEAM    = {header.get('BEAM')}")
    print(f"  BMAJ    = {header.get('BMAJ')}")
    print(f"  BMIN    = {header.get('BMIN')}")

    if data is None:
        print("  No data in this HDU.")
        return

    if not isinstance(data, np.ndarray):
        print("  Data exists but is not a normal numpy image array.")
        return

    if not np.issubdtype(data.dtype, np.number):
        print(f"  Data exists but is not numeric image data. dtype = {data.dtype}")
        return

    arr = np.asarray(data, dtype=float)
    finite = np.isfinite(arr)

    print(f"  DATA SHAPE = {arr.shape}")
    print(f"  DATA DTYPE = {arr.dtype}")
    print(f"  FINITE PIXELS = {finite.sum()} / {arr.size}")

    if finite.sum() > 0:
        print(f"  MIN    = {np.nanmin(arr):.6e}")
        print(f"  MAX    = {np.nanmax(arr):.6e}")
        print(f"  MEAN   = {np.nanmean(arr):.6e}")
        print(f"  MEDIAN = {np.nanmedian(arr):.6e}")


print("\n==============================")
print("RAW HERSCHEL FITS FILES")
print("==============================")

for file in fits_files:
    print(f"\nFILE: {file.name}")
    print("-" * 80)

    with fits.open(file) as hdul:
        hdul.info()

        for i, hdu in enumerate(hdul):
            inspect_hdu(hdu, i)


print("\n==============================")
print("KERNEL FITS FILES")
print("==============================")

for file in kernel_files:
    print(f"\nKERNEL: {file.name}")
    print("-" * 80)

    with fits.open(file) as hdul:
        hdul.info()

        for i, hdu in enumerate(hdul):
            inspect_hdu(hdu, i)