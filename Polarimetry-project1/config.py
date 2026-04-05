from pathlib import Path

# Base paths
BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent

# Input directories
BIAS_DIR = PROJECT_ROOT / "bias"
SCIENCE_DIR = PROJECT_ROOT / "HD212311" / "HD212311"

# Output directory
OUTPUT_DIR = BASE_DIR / "output"

# Subdirectories
MASTER_BIAS_DIR = OUTPUT_DIR / "master_bias"
BIAS_SUB_DIR = OUTPUT_DIR / "bias_subtracted"
STACK_DIR = OUTPUT_DIR / "stacks"
LOG_DIR = OUTPUT_DIR / "logs"
PLOT_DIR = OUTPUT_DIR / "plots"
PHOT_DIR = OUTPUT_DIR / "photometry"
POL_DIR = OUTPUT_DIR / "polarimetry"

OUTPUT_SUBDIRS = [
    OUTPUT_DIR,
    MASTER_BIAS_DIR,
    BIAS_SUB_DIR,
    STACK_DIR,
    LOG_DIR,
    PLOT_DIR,
    PHOT_DIR,
    POL_DIR,
]

# HWP mapping
GROUP_TO_HWP = {
    1: 67.5,
    2: 45.0,
    3: 22.5,
    4: 0.0,
}

# Detection parameters
DAO_FWHM = 4.0
DAO_THRESHOLD_SIGMA = 5.0

# Photometry parameters
APERTURE_RADII = [3, 4, 5, 6, 8]
ANNULUS_R_IN = 10
ANNULUS_R_OUT = 15

# Beam pairing (initial guess, adjust later)
DX_BEAM = 31.0
DY_BEAM = 5.0
PAIR_TOL = 4.0

# Alignment star selection
ALIGNMENT_GROUP_ID = 1
ALIGNMENT_X = 102
ALIGNMENT_Y = 133
ALIGNMENT_RADIUS = 6
ALIGNMENT_CUTOUT_HALF_SIZE = 12