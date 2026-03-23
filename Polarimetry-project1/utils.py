import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from astropy.io import fits
from astropy.visualization import simple_norm


def ensure_directories(dir_list):
    for d in dir_list:
        d.mkdir(parents=True, exist_ok=True)


def clean_header(header):
    """
    Return a cleaned FITS header by keeping only cards that astropy can write.
    Corrupted cards are dropped.
    """
    clean = fits.Header()

    for card in header.cards:
        try:
            clean.append(card)
        except Exception:
            continue

    return clean


def load_fits_data(path):
    """
    Load FITS image safely and squeeze extra dimensions.
    Returns data and a cleaned header.
    """
    with fits.open(path, ignore_missing_simple=True) as hdul:
        data = hdul[0].data
        raw_header = hdul[0].header.copy()

    if data is None:
        raise ValueError(f"No image data found in {path}")

    data = np.squeeze(data).astype(float)
    header = clean_header(raw_header)

    return data, header


def save_fits(path, data, header=None):
    """
    Save FITS safely.

    If header is bad or missing, write without it.
    """
    data = np.asarray(data, dtype=np.float32)

    if header is None:
        fits.writeto(path, data, overwrite=True, output_verify="silentfix")
        return

    try:
        header = clean_header(header)
        fits.writeto(path, data, header=header, overwrite=True, output_verify="silentfix")
    except Exception:
        fits.writeto(path, data, overwrite=True, output_verify="silentfix")


def show_and_save_image(data, outpath=None, title="", percentile=99.5, cmap="gray"):
    norm = simple_norm(data, stretch="linear", percent=percentile)
    plt.figure(figsize=(7, 7))
    plt.imshow(data, origin="lower", cmap=cmap, norm=norm)
    plt.colorbar(label="Counts")
    plt.xlabel("X pixel")
    plt.ylabel("Y pixel")
    plt.title(title)
    plt.tight_layout()

    if outpath is not None:
        plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()


def parse_science_filename(path):
    """
    Parse names like hd21p11.fits
    Returns (group_id, exp_id)
    """
    m = re.search(r"p(\d)(\d)\.fits$", path.name.lower())
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def build_science_dataframe(science_files):
    rows = []
    for f in science_files:
        group_id, exp_id = parse_science_filename(f)
        rows.append((f.name, group_id, exp_id))

    rows = sorted(
        rows,
        key=lambda x: (
            999 if x[1] is None else x[1],
            999 if x[2] is None else x[2],
            x[0]
        )
    )

    return pd.DataFrame(rows, columns=["file", "group_id", "exp_id"])


def plot_marked_star(data, x, y, radius=12, title="", outpath=None, cmap="gray"):
    """
    Plot full image with a circle marking the selected star.
    """
    norm = simple_norm(data, stretch="linear", percent=99.5)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(data, origin="lower", cmap=cmap, norm=norm)

    circ = Circle((x, y), radius, edgecolor="red", facecolor="none", linewidth=2)
    ax.add_patch(circ)

    ax.plot(x, y, marker="+", color="yellow", markersize=10, mew=2)

    ax.set_xlabel("X pixel")
    ax.set_ylabel("Y pixel")
    ax.set_title(title)

    if outpath is not None:
        plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()


def extract_cutout(data, x, y, half_size=15):
    """
    Return a square cutout centered near (x, y).
    """
    x = int(round(x))
    y = int(round(y))

    x1 = max(0, x - half_size)
    x2 = min(data.shape[1], x + half_size + 1)
    y1 = max(0, y - half_size)
    y2 = min(data.shape[0], y + half_size + 1)

    return data[y1:y2, x1:x2], x1, y1


def plot_cutout_with_circle(data, x, y, radius=12, half_size=15, title="", outpath=None, cmap="gray"):
    """
    Plot zoomed cutout around selected star with circle and center mark.
    """
    cutout, x1, y1 = extract_cutout(data, x, y, half_size=half_size)

    x_local = x - x1
    y_local = y - y1

    norm = simple_norm(cutout, stretch="linear", percent=99.5)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(cutout, origin="lower", cmap=cmap, norm=norm)

    circ = Circle((x_local, y_local), radius, edgecolor="red", facecolor="none", linewidth=2)
    ax.add_patch(circ)

    ax.plot(x_local, y_local, marker="+", color="yellow", markersize=10, mew=2)

    ax.set_xlabel("X pixel in cutout")
    ax.set_ylabel("Y pixel in cutout")
    ax.set_title(title)

    if outpath is not None:
        plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()