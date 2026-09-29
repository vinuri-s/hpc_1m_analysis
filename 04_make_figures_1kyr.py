"""
Analysis figures for the 1 m native-resolution, 1000yr/DT=10yr sensitivity
run (results_1kyr_dt10/) -- this is now the PRODUCTION run for the paper
(replacing the 1 Myr/DT=5000yr run, which showed a numerical instability
in SPACE at K6/K7 traced to the coarser timestep; see 03_run_one_1m_1kyr.py
header). Run this AFTER all 14 configs have finished and written
results_1kyr_dt10/{model}_{kname}.npz + _log.json.

Produces the same three figure types as the original 1 Myr analysis
(04_make_figures_1m.py), relabelled for the 1000yr/10yr run so they can
be dropped into the paper in place of the old ones:
  - fig_sensitivity_timeseries_1kyr.png   (2x4 grid: relief/slope/erosion/mean-z)
  - fig_diffmaps_by_K_1kyr.png            (2x7 grid: elevation change per K)
  - fig_diffmaps_structural_1kyr.png      (1x7 grid: SPACE-SPIM difference)

Usage:
    python3 04_make_figures_1kyr.py
"""
import io
import json
import pickle
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import SymLogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
RES_DIR = os.path.join(HERE, "results_1kyr_dt10")


def robust_read_bytes(path, retries=60, delay=1.0):
    """This OneDrive-backed mount intermittently raises EDEADLK (Resource
    deadlock avoided) on open()/read() -- including partway through a
    read, so np.load()'s lazy per-array reads from an on-disk .npz can
    fail even after the file "opened" successfully. Work around this by
    retrying the *entire* open+read(all bytes)+close as one atomic unit,
    so every caller downstream (json.loads, np.load(BytesIO(...))) then
    operates purely on an in-memory buffer with no further filesystem
    access for that file.
    """
    last_err = None
    for attempt in range(retries):
        try:
            with open(path, "rb") as f:
                return f.read()
        except OSError as e:
            last_err = e
            print(f"    [retry {attempt+1}/{retries}] {os.path.basename(path)}: {e}", flush=True)
            time.sleep(delay)
    raise last_err


def robust_load_npz(path):
    return np.load(io.BytesIO(robust_read_bytes(path)))


base = pickle.loads(robust_read_bytes(os.path.join(HERE, "grid_base_1m.pkl")))
ny, nx = base["shape"]
res = base["res"]
valid = base["valid_flat"].reshape(ny, nx)
z0 = base["z0"].reshape(ny, nx)
core = (base["status"] == 0).reshape(ny, nx)  # BC_NODE_IS_CORE == 0

# UCS range 3.55-217.02 MPa (99th percentile) -- identical K values to the
# original 1 Myr sweep; see 03_run_one_1m_1kyr.py header for the
# derivation.
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
# Okabe-Ito colorblind-safe qualitative palette (Okabe & Ito, 2008), one
# clearly distinguishable colour per K value (K1 = most resistant, blue,
# through K7 = most erodible, reddish purple), replacing the previous
# viridis_r sequential colormap whose middle values (K3-K5) were hard to
# tell apart at a glance.
OKABE_ITO_7 = ["#0072B2", "#56B4E9", "#009E73", "#F0E442",
               "#E69F00", "#D55E00", "#CC79A7"]
COLORS = {k: OKABE_ITO_7[i] for i, k in enumerate(KEYS_ORDERED)}
LABELS = {k: rf"$K_{{br}}={K_VALUES[k]:.2e}$" for k in KEYS_ORDERED}

DT = 10.0  # yr

# ---------------------------------------------------------------------
# Load logs (cheap, kept in memory for the whole script) and snapshots
# (expensive: 14 configs x ~11 snapshots x 7.2M nodes would need >13 GB
# held simultaneously, which does not fit in this machine's RAM). Instead,
# each config's .npz is opened once, its relief time series is computed
# immediately, and only the final-state array (~58 MB) is cached before
# the full snapshot dict is dropped -- peak memory stays under ~2 GB.
# ---------------------------------------------------------------------
# Local (non-OneDrive) cache so that repeated invocations of this script
# -- needed because the OneDrive-backed RES_DIR mount intermittently
# deadlocks (EDEADLK) on individual files -- can resume from wherever the
# last run left off instead of re-reading every already-processed config
# from the flaky mount each time.
import tempfile
CACHE_DIR = os.path.join(tempfile.gettempdir(), "fig_1kyr_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

logs = {}
finals = {}   # finals[(model, k)] -> final-timestep elevation array (ny, nx)
relief_ts = {}
STEPS = None
for model in ["spim", "space"]:
    logs[model] = {}
    relief_ts[model] = {}
    for k in KEYS_ORDERED:
        cache_path = os.path.join(CACHE_DIR, f"{model}_{k}.pkl")
        if os.path.exists(cache_path):
            with open(cache_path, "rb") as f:
                cached = pickle.load(f)
            logs[model][k] = cached["log"]
            relief_ts[model][k] = cached["relief_ts"]
            finals[(model, k)] = cached["final"]
            if STEPS is None:
                STEPS = cached["steps"]
            print(f"  [cache] loaded {model}_{k}", flush=True)
            continue
        logs[model][k] = json.loads(robust_read_bytes(f"{RES_DIR}/{model}_{k}_log.json"))
        d = robust_load_npz(f"{RES_DIR}/{model}_{k}.npz")
        steps = sorted(int(s.split("_")[1]) for s in d.files if s.startswith("snap_"))
        if STEPS is None:
            STEPS = steps
        vals = []
        for s in steps:
            zc = d[f"snap_{s}"].reshape(ny, nx)[core]
            vals.append(np.percentile(zc, 90) - np.percentile(zc, 10))
        relief_ts[model][k] = np.array(vals)
        finals[(model, k)] = d[f"snap_{steps[-1]}"].reshape(ny, nx).copy()
        d.close()
        with open(cache_path, "wb") as f:
            pickle.dump({"log": logs[model][k], "relief_ts": relief_ts[model][k],
                         "final": finals[(model, k)], "steps": steps}, f)
        print(f"  processed {model}_{k}", flush=True)

FINAL = STEPS[-1]
t_snap = np.array(STEPS) * DT  # years (DT=10 yr, total 1000 yr)

def add_centered_row(fig, gs, row, n_items, ncols):
    """Add n_items subplots into gridspec row `row` (which is divided into
    ncols*2 fine columns, 2 per item), centering the row when n_items <
    ncols instead of packing items flush to the left. This replaces the
    previous behaviour of leaving a ragged blank panel on the right of an
    incomplete row.
    """
    left_pad = ncols - n_items
    row_axes = []
    for i in range(n_items):
        start = left_pad + 2 * i
        ax = fig.add_subplot(gs[row, start:start + 2])
        row_axes.append(ax)
    return row_axes


def add_vertical_block_label(fig, top_ax, bottom_ax, text, x):
    """Add a bold, rotated label vertically centred on the block spanned by
    top_ax (top row of the block) through bottom_ax (bottom row), at a
    fixed figure-fraction x position -- used for the SPIM/SPACE row
    labels, which previously sat at the top of their block instead of
    being centred on it. `x` must be chosen outside the margin already
    occupied by each panel's own y-axis label/ticks, otherwise the two
    overlap.
    """
    bbox_top = top_ax.get_position()
    bbox_bot = bottom_ax.get_position()
    y_center = (bbox_top.y1 + bbox_bot.y0) / 2
    fig.text(x, y_center, text, rotation=90,
              va="center", ha="center", fontsize=12, fontweight="bold")


# A4 width is 8.27 in; keep every figure at or under that so it drops
# straight into the journal page without downscaling.
A4_WIDTH_IN = 8.27

# ---------------------------------------------------------------------
# (1) Sensitivity time series. Previously a 2x4 grid (one row per model,
# 4 metrics across), which put four panels in a single row. Restructured
# as two stacked 2x2 blocks (one per model, 2 metrics per row) so each
# row now holds only 2 panels -- consistent with the block layout used
# below for the diff-map figures, and with each panel noticeably larger.
# ---------------------------------------------------------------------
NCOLS_TS = 2
metrics = [
    ("relief", "P90-P10 relief (m)"),
    ("mean_slope", "Mean slope (m/m)"),
    ("mean_erosion_rate", "Mean erosion rate (m/yr)\n(+ = net erosion, - = net aggradation)"),
    ("mean_elev_core", "Mean core elevation (m)"),
]
# Slightly wider than A4_WIDTH_IN here specifically: this is a source
# raster that LaTeX will scale down to fit the journal column width
# anyway (aspect ratio and relative panel/gap proportions are what
# actually matter, not the literal inch count), and the extra width lets
# the two columns keep a generous gap (see wspace below) without
# shrinking the panels themselves to make room for it.
fig, axes = plt.subplots(4, NCOLS_TS, figsize=(A4_WIDTH_IN + 1.0, 13.0), sharex=True)
model_block_axes = {}
for m_idx, model in enumerate(["spim", "space"]):
    block_row0 = m_idx * 2
    block_axes = []
    for metric_idx, (metric, ylabel) in enumerate(metrics):
        r_in_block, col = divmod(metric_idx, NCOLS_TS)
        row = block_row0 + r_in_block
        ax = axes[row, col]
        block_axes.append(ax)
        for k in KEYS_ORDERED:
            if metric == "relief":
                t, y = t_snap, relief_ts[model][k]
            else:
                t = np.array(logs[model][k]["t"])
                y = np.array(logs[model][k][metric])
            ax.plot(t, y, color=COLORS[k], label=LABELS[k], lw=1.6,
                     marker="o" if metric == "relief" else None, ms=3)
        if metric == "mean_erosion_rate":
            ax.axhline(0, color="grey", lw=0.5, ls=":")
        if row == 3:
            ax.set_xlabel("Model time (yr)")
        ax.set_ylabel(ylabel, fontsize=8.5)
        ax.set_title(chr(97 + row * NCOLS_TS + col) + ")", loc="left",
                      fontweight="bold", fontsize=10)
    model_block_axes[model] = block_axes

plt.suptitle("1 m native-resolution sensitivity sweep, 1000 yr / DT=10 yr\n"
             "production run: StreamPowerEroder (top) vs SpaceLargeScaleEroder\n"
             "with $K_{sed}=10\\times K_{br}$ (bottom)",
             fontsize=11)
plt.tight_layout(rect=[0.0, 0.0, 0.82, 0.97])
# Reserve a clear left-hand strip for the SPIM/SPACE block labels, outside
# the margin tight_layout already gave the per-panel y-axis labels, so the
# two don't overlap, and leave a little breathing room on both sides of
# the label text itself so it doesn't sit flush against the figure edge.
# Also widen the gap between the two columns -- the narrower panel width
# from reserving space on the right for the colour-key colorbar left too
# little room between e.g. panel (e)'s tick labels and panel (f)'s y-axis
# label, so they were running into each other.
fig.subplots_adjust(left=fig.subplotpars.left + 0.13, wspace=0.45)
fig.canvas.draw()
for model in ["spim", "space"]:
    block_axes = model_block_axes[model]
    label = "SPIM\n(StreamPowerEroder)" if model == "spim" else "SPACE\n(SpaceLargeScaleEroder)"
    add_vertical_block_label(fig, block_axes[0], block_axes[-1], label, x=0.045)

# Discrete K1-K7 colour key, styled like the colorbars on the diff-map
# figures below, instead of a bottom-of-figure line legend -- keeps the
# colour -> K mapping consistent in presentation across all figures.
from matplotlib.colors import ListedColormap, BoundaryNorm
ts_cmap = ListedColormap(OKABE_ITO_7)
ts_norm = BoundaryNorm(np.arange(len(KEYS_ORDERED) + 1), ts_cmap.N)
ts_sm = plt.cm.ScalarMappable(cmap=ts_cmap, norm=ts_norm)
ts_cbar_ax = fig.add_axes([0.86, 0.30, 0.03, 0.40])
ts_cbar = fig.colorbar(ts_sm, cax=ts_cbar_ax,
                        ticks=np.arange(len(KEYS_ORDERED)) + 0.5)
ts_cbar.ax.set_yticklabels([LABELS[k] for k in KEYS_ORDERED], fontsize=8)
ts_cbar.set_label(r"$K_{br}$ (yr$^{-1}$)", fontsize=9)

plt.savefig(os.path.join(HERE, "fig_sensitivity_timeseries_1kyr.png"), dpi=260,
            bbox_inches="tight")
print("Saved fig_sensitivity_timeseries_1kyr.png")

# ---------------------------------------------------------------------
# (2) Within-model diff maps. Two stacked blocks (one per model), each
# holding its 7 K panels as a full row of 4 plus a second row of 3. The
# row of 3 is now centred under the row of 4 (via a fine gridspec, 2
# sub-columns per panel) rather than left-packed with a ragged blank
# panel, and the SPIM/SPACE row labels are centred vertically on their
# 2-row block instead of sitting at its top.
# ---------------------------------------------------------------------
NCOLS = 4
FCOLS = NCOLS * 2
fig = plt.figure(figsize=(A4_WIDTH_IN, 8.7))
gs = fig.add_gridspec(4, FCOLS, hspace=0.12, wspace=0.05)
extent = [0, nx * res, 0, ny * res]

dz_all = {}
for model in ["spim", "space"]:
    for k in KEYS_ORDERED:
        dz_all[(model, k)] = finals[(model, k)] - z0

vmax = np.nanmax([np.nanpercentile(np.abs(dz_all[key][valid]), 99) for key in dz_all])
vmax = max(vmax, 1e-6)
linthresh = max(vmax * 0.01, 1e-6)  # auto-scaled, unlike the fixed 2.0 m
                                     # used for the 1 Myr run's much larger dz
norm = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)

all_axes = []
model_block_axes = {}
for m_idx, model in enumerate(["spim", "space"]):
    block_row0 = m_idx * 2
    row0_keys = KEYS_ORDERED[:NCOLS]
    row1_keys = KEYS_ORDERED[NCOLS:]
    row_axes = (add_centered_row(fig, gs, block_row0, len(row0_keys), NCOLS)
                + add_centered_row(fig, gs, block_row0 + 1, len(row1_keys), NCOLS))
    for ax, k in zip(row_axes, KEYS_ORDERED):
        dz = np.ma.masked_where(~valid, dz_all[(model, k)])
        im = ax.imshow(dz, cmap="RdBu_r", norm=norm, origin="upper", extent=extent)
        ax.set_title(rf"$K_{{br}}$={K_VALUES[k]:.2e}", fontsize=9)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
    model_block_axes[model] = row_axes
    all_axes.extend(row_axes)

plt.suptitle("1 m native-resolution elevation-change maps after 1000 yr "
             "(uplift-only baseline = +0.008 m everywhere)\n"
             "top block: SPIM; bottom block: SPACE", fontsize=12)
# Panels are added via GridSpec with per-row column spans that vary (to
# centre the row of 3), which tight_layout cannot reason about -- set
# margins directly instead. The colorbar gets its own explicitly placed
# axes (rather than fig.colorbar(..., ax=all_axes), which steals space
# from the panel axes at the time it's called and collided with the
# panels once the margins above were applied afterwards).
fig.subplots_adjust(left=0.11, right=0.78, top=0.86, bottom=0.02)
cbar_ax = fig.add_axes([0.83, 0.15, 0.03, 0.55])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(r"$\Delta z$ over 1000 yr (m, symmetric-log scale); "
                "blue = net erosion, red = net aggradation", fontsize=9)
fig.canvas.draw()
for model in ["spim", "space"]:
    block_axes = model_block_axes[model]
    add_vertical_block_label(fig, block_axes[0], block_axes[-1], model.upper(), x=0.04)
plt.savefig(os.path.join(HERE, "fig_diffmaps_by_K_1kyr.png"), dpi=260, bbox_inches="tight")
print("Saved fig_diffmaps_by_K_1kyr.png")

# ---------------------------------------------------------------------
# (3) Structural (SPACE - SPIM) diff maps. A row of 4 plus a row of 3,
# with the row of 3 centred under the row of 4 (same fine-gridspec
# approach as the diff-by-K figure above) instead of left-packed with a
# ragged blank panel.
# ---------------------------------------------------------------------
fig = plt.figure(figsize=(A4_WIDTH_IN, 4.6))
gs = fig.add_gridspec(2, FCOLS, hspace=0.15, wspace=0.05)
struct_diff = {k: finals[("space", k)] - finals[("spim", k)] for k in KEYS_ORDERED}

vmax2 = np.nanmax([np.nanpercentile(np.abs(struct_diff[k][valid]), 99) for k in KEYS_ORDERED])
vmax2 = max(vmax2, 1e-6)
linthresh2 = max(vmax2 * 0.01, 1e-6)
norm2 = SymLogNorm(linthresh=linthresh2, vmin=-vmax2, vmax=vmax2, base=10)

row0_keys = KEYS_ORDERED[:NCOLS]
row1_keys = KEYS_ORDERED[NCOLS:]
struct_axes = (add_centered_row(fig, gs, 0, len(row0_keys), NCOLS)
               + add_centered_row(fig, gs, 1, len(row1_keys), NCOLS))
for ax, k in zip(struct_axes, KEYS_ORDERED):
    dd = np.ma.masked_where(~valid, struct_diff[k])
    im = ax.imshow(dd, cmap="PuOr_r", norm=norm2, origin="upper", extent=extent)
    ax.set_title(rf"$K_{{br}}$={K_VALUES[k]:.2e}", fontsize=9)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])

plt.suptitle("1 m native-resolution structural sensitivity: SPACE vs SPIM "
             "at matched bedrock erodibility ($K_{sed}=10\\times K_{br}$)", fontsize=11)
# Same reason as the diff-by-K figure above: variable per-row column
# spans aren't compatible with tight_layout, so set margins directly, and
# give the colorbar its own explicitly placed axes rather than letting
# fig.colorbar(..., ax=struct_axes) steal space from the panels at call
# time (which then collided with the panels once these margins were
# applied afterwards).
fig.subplots_adjust(left=0.05, right=0.78, top=0.80, bottom=0.03)
cbar_ax = fig.add_axes([0.83, 0.20, 0.03, 0.55])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(r"$z_{SPACE}-z_{SPIM}$ at 1000 yr (m, symmetric-log scale)", fontsize=9)
plt.savefig(os.path.join(HERE, "fig_diffmaps_structural_1kyr.png"), dpi=260, bbox_inches="tight")
print("Saved fig_diffmaps_structural_1kyr.png")

print("\n=== SUMMARY (final state, t=1000 yr) ===")
for model in ["spim", "space"]:
    for k in KEYS_ORDERED:
        L = logs[model][k]
        print(f"{model:6s} {k:14s} K={K_VALUES[k]:.3e}  "
              f"relief={relief_ts[model][k][-1]:.4f} m  "
              f"mean_slope={L['mean_slope'][-1]:.6f}  "
              f"mean_eros={L['mean_erosion_rate'][-1]:.3e} m/yr")
for k in KEYS_ORDERED:
    dd = struct_diff[k][valid]
    print(f"SPACE-SPIM ({k}): mean={np.nanmean(dd):.5f} m, "
          f"min={np.nanmin(dd):.5f}, max={np.nanmax(dd):.5f}")
