from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PLOT_DIR = ROOT / "scripts" / "plots"
PLOT_DIR.mkdir(exist_ok=True)

file = ROOT / "ngc1866_ms" / "J_AJ_118_2839_stars.csv"
df = pd.read_csv(file)

df = df.dropna(subset=["Xpos", "Ypos", "Vmag", "B-V", "e_Vmag", "e_B-V"]).copy()

# Basic photometric quality cut
df = df[
    (df["Vmag"] > 10) & (df["Vmag"] < 24) &
    (df["B-V"] > -0.5) & (df["B-V"] < 2.0) &
    (df["e_Vmag"] < 0.15) &
    (df["e_B-V"] < 0.20)
].copy()

# Estimate cluster center from median position
x0 = df["Xpos"].median()
y0 = df["Ypos"].median()

df["r_pix"] = np.sqrt((df["Xpos"] - x0)**2 + (df["Ypos"] - y0)**2)

print("Estimated center:")
print("X0 =", x0, "Y0 =", y0)

print("\nRadius statistics:")
print(df["r_pix"].describe())

# Choose core and field regions
# You may adjust these after seeing the spatial plot.
core = df[df["r_pix"] < 350].copy()
field = df[df["r_pix"] > 650].copy()

print("\nTotal stars:", len(df))
print("Core stars:", len(core))
print("Outer field stars:", len(field))

# -----------------------------
# 1. Spatial distribution
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(df["Xpos"], df["Ypos"], s=2, alpha=0.20, label="All stars")
plt.scatter(core["Xpos"], core["Ypos"], s=3, alpha=0.55, label="Cluster core")
plt.xlabel("X position [pixel]")
plt.ylabel("Y position [pixel]")
plt.title("NGC 1866: Spatial Selection")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "ngc1866_spatial_selection.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 2. CMD core vs field
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(field["B-V"], field["Vmag"], s=3, alpha=0.15, label="Outer field")
plt.scatter(core["B-V"], core["Vmag"], s=4, alpha=0.45, label="Cluster core")

plt.gca().invert_yaxis()
plt.xlabel("B - V")
plt.ylabel("V")
plt.title("NGC 1866 CMD: Core vs Field")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "ngc1866_core_field_cmd.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 3. Zoom on main sequence
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(field["B-V"], field["Vmag"], s=3, alpha=0.10, label="Outer field")
plt.scatter(core["B-V"], core["Vmag"], s=4, alpha=0.55, label="Cluster core")

plt.xlim(-0.3, 0.8)
plt.ylim(23.5, 14.0)
plt.xlabel("B - V")
plt.ylabel("V")
plt.title("NGC 1866 Main Sequence Region")
plt.legend()
plt.grid(True)
plt.savefig(PLOT_DIR / "ngc1866_main_sequence_zoom.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------
# 4. Compare expected LMC distance modulus
# -----------------------------
mu_lmc = 18.49
E_BV = 0.08
A_V = 3.1 * E_BV

print("\nExpected LMC comparison:")
print("mu_LMC =", mu_lmc)
print("Assumed E(B-V) =", E_BV)
print("A_V =", A_V)
print("Apparent shift in V = mu + A_V =", mu_lmc + A_V)

with open(ROOT / "ngc1866_ms" / "ngc1866_spatial_cmd_notes.txt", "w") as f:
    f.write("NGC 1866 spatial CMD analysis\n")
    f.write(f"Total stars = {len(df)}\n")
    f.write(f"Core stars = {len(core)}\n")
    f.write(f"Field stars = {len(field)}\n")
    f.write(f"Estimated center: X0={x0:.2f}, Y0={y0:.2f}\n")
    f.write("No VPD selection possible because this catalogue lacks proper motions.\n")
    f.write("Use Xpos/Ypos spatial selection instead.\n")
    f.write(f"Expected LMC distance modulus = {mu_lmc:.2f}\n")
    f.write(f"Assumed E(B-V) = {E_BV:.2f}, A_V = {A_V:.2f}\n")

print("\nSaved plots in:", PLOT_DIR)