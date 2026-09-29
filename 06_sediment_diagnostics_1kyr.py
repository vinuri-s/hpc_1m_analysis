"""
1000yr/DT=10yr version of 06_sediment_diagnostics.py: SPACE soil-depth
maps, sediment-flux/soil-depth time series, and a mass-balance sanity
check, for the results_1kyr_dt10/ production run.

Usage:
    python3 06_sediment_diagnostics_1kyr.py
"""
import os
import json
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
# Okabe-Ito colorblind-safe qualitative palette (Okabe & Ito, 2008), same
# K1-K7 assignment used across all figures in the paper (see
# 04_make_figures_1kyr.py), replacing viridis_r so K colours are
# consistent between the main-text and supplementary figures.
OKABE_ITO_7 = ["#0072B2", "#56B4E9", "#009E73", "#F0E442",
               "#E69F00", "#D55E00", "#CC79A7"]
COLORS = {k: OKABE_ITO_7[i] for i, k in enumerate(KEYS_ORDERED)}

logs = {}
soil_final = {}
for k in KEYS_ORDERED:
    with open(f"{RES_DIR}/space_{k}_log.json") as f:
        logs[k] = json.load(f)
    d = np.load(f"{RES_DIR}/space_{k}.npz")
    soil_steps = sorted(int(s.split("_")[1]) for s in d.files if s.startswith("soil_"))
    soil_final[k] = d[f"soil_{soil_steps[-1]}"].reshape(ny, nx)

fig, axes = plt.subplots(1, 7, figsize=(24, 4.2))
extent = [0, nx * res, 0, ny * res]
vmax = np.nanmax([np.nanpercentile(soil_final[k][valid], 99) for k in KEYS_ORDERED])

for ax, k in zip(axes, KEYS_ORDERED):
    s = np.ma.masked_where(~valid, soil_final[k])
    im = ax.imshow(s, cmap="YlOrBr", vmin=0, vmax=vmax, origin="upper", extent=extent)
    ax.set_title(rf"$K_{{br}}$={K_VALUES[k]:.2e}", fontsize=9)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])

cbar = fig.colorbar(im, ax=axes, shrink=0.75, pad=0.01)
cbar.set_label("Soil (regolith) depth at 1000 yr (m)")
plt.suptitle("1 m, 1000yr/DT=10yr: SPACE sediment cover, final soil depth across the basalt erodibility range", fontsize=12)
plt.savefig(f"{HERE}/fig_soil_depth_maps_1kyr.png", dpi=220, bbox_inches="tight")
print("Saved fig_soil_depth_maps_1kyr.png")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
for k in KEYS_ORDERED:
    t = np.array(logs[k]["t"])  # years
    axes[0].plot(t, logs[k]["mean_soil_depth"], color=COLORS[k], lw=1.6,
                 label=rf"$K_{{br}}={K_VALUES[k]:.2e}$")
    axes[1].plot(t, logs[k]["sediment_flux_outlet"], color=COLORS[k], lw=1.6)
axes[0].axhline(1.0, color="grey", lw=0.5, ls=":")
axes[0].set_xlabel("Model time (yr)"); axes[0].set_ylabel("Mean soil depth (m)")
axes[1].set_xlabel("Model time (yr)"); axes[1].set_ylabel(r"Sediment flux at outlet (m$^3$/yr)")
axes[0].legend(fontsize=7, loc="best")
plt.tight_layout()
plt.savefig(f"{HERE}/fig_sediment_flux_soil_ts_1kyr.png", dpi=230, bbox_inches="tight")
print("Saved fig_sediment_flux_soil_ts_1kyr.png")

print("\n=== Mass-balance sanity check (SPACE, 1 m, 1000yr/DT=10yr) ===")
cell_area = res * res
for k in KEYS_ORDERED:
    d = np.load(f"{RES_DIR}/space_{k}.npz")
    steps = sorted(int(s.split("_")[1]) for s in d.files if s.startswith("snap_"))
    z0 = d[f"snap_{steps[0]}"].reshape(ny, nx)
    zf = d[f"snap_{steps[-1]}"].reshape(ny, nx)
    t = np.array(logs[k]["t"])
    uplift_total = 8.0e-6 * (t[-1] - t[0])
    dz = (zf - z0)[valid] - uplift_total
    vol_change = np.sum(dz) * cell_area
    qs = np.array(logs[k]["sediment_flux_outlet"])
    vol_exported = np.trapezoid(qs, t)
    print(f"{k}: vol change (excl. uplift) = {vol_change:,.3f} m^3   "
          f"cumulative exported = {vol_exported:,.3f} m^3   "
          f"ratio = {vol_exported/abs(vol_change) if vol_change != 0 else float('nan'):.2f}")
