"""
1000yr/DT=10yr version of 08_final_state_maps.py: hillshaded
final-elevation maps for each K value and model, for the
results_1kyr_dt10/ production run.

Note: because total imposed change over 1000 yr is only mm-cm scale
against a real DEM with hundreds of metres of relief, these seven panels
per model will look visually almost identical to the t=0 topography and
to each other -- that is the expected, correct result at this duration,
not a bug. The quantitative differences are in fig_diffmaps_by_K_1kyr.png
and fig_sensitivity_timeseries_1kyr.png instead.

Usage:
    python3 08_final_state_maps_1kyr.py
"""
import os
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

HERE = os.path.dirname(os.path.abspath(__file__))
RES_DIR = os.path.join(HERE, "results_1kyr_dt10")

with open(os.path.join(HERE, "grid_base_1m.pkl"), "rb") as f:
    base = pickle.load(f)
ny, nx = base["shape"]
res = base["res"]
valid = base["valid_flat"].reshape(ny, nx)

K_VALUES = {
    "K1_3.02e-09": 3.024e-09,
    "K2_1.19e-08": 1.191e-08,
    "K3_4.69e-08": 4.693e-08,
    "K4_1.85e-07": 1.849e-07,
    "K5_7.28e-07": 7.283e-07,
    "K6_2.87e-06": 2.869e-06,
    "K7_1.13e-05": 1.130e-05,
}
KEYS_ORDERED = list(K_VALUES.keys())
extent = [0, nx * res, 0, ny * res]
ls = LightSource(azdeg=315, altdeg=45)

fig, axes = plt.subplots(2, 7, figsize=(24, 7.2))
zfinal = {}
for model in ["spim", "space"]:
    for k in KEYS_ORDERED:
        d = np.load(f"{RES_DIR}/{model}_{k}.npz")
        steps = sorted(int(s.split("_")[1]) for s in d.files if s.startswith("snap_"))
        zfinal[(model, k)] = d[f"snap_{steps[-1]}"].reshape(ny, nx)

vmin = np.nanmin([np.nanmin(np.where(valid, zfinal[key], np.nan)) for key in zfinal])
vmax = np.nanmax([np.nanmax(np.where(valid, zfinal[key], np.nan)) for key in zfinal])

for row, model in enumerate(["spim", "space"]):
    for col, k in enumerate(KEYS_ORDERED):
        ax = axes[row, col]
        z = zfinal[(model, k)]
        z_fill = np.where(valid, z, np.nanmin(z[valid]))
        rgb = ls.shade(z_fill, cmap=plt.cm.terrain, vmin=vmin, vmax=vmax,
                        blend_mode="soft", vert_exag=2.0)
        ax.imshow(rgb, origin="upper", extent=extent)
        ax.imshow(np.ma.masked_where(valid, np.ones_like(z)), cmap="Greys",
                  vmin=0, vmax=1, origin="upper", extent=extent)
        if row == 0:
            ax.set_title(rf"$K_{{br}}$={K_VALUES[k]:.2e}", fontsize=9)
        if col == 0:
            ax.set_ylabel("SPIM" if row == 0 else "SPACE", fontsize=11, fontweight="bold")
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])

sm = plt.cm.ScalarMappable(cmap="terrain", norm=plt.Normalize(vmin=vmin, vmax=vmax))
cbar = fig.colorbar(sm, ax=axes, shrink=0.7, pad=0.01)
cbar.set_label("Elevation at t=1000 yr (m)")
plt.suptitle("1 m: final landscape morphology (hillshaded elevation) after 1000 yr / DT=10yr", fontsize=13)
plt.savefig(f"{HERE}/fig_final_state_maps_1kyr.png", dpi=220, bbox_inches="tight")
print("Saved fig_final_state_maps_1kyr.png")
