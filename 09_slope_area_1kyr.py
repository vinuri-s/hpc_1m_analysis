"""
1000yr/DT=10yr version of 09_slope_area.py: binned slope-area regression
with fitted concavity (theta) and steepness index (ks) per K value/model,
final-state (t=1000 yr) topography, for the results_1kyr_dt10/ production
run.

Usage:
    python3 09_slope_area_1kyr.py
"""
import os
import ctypes
from ctypes import wintypes
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, DepressionFinderAndRouter

HERE = os.path.dirname(os.path.abspath(__file__))
RES_DIR = os.path.join(HERE, "results_1kyr_dt10")

STACK_SIZE_BYTES = 1024 * 1024 * 1024  # 1 GiB

kernel32 = ctypes.windll.kernel32
THREAD_ROUTINE = ctypes.WINFUNCTYPE(wintypes.DWORD, wintypes.LPVOID)


def run_in_big_stack(func, *args, **kwargs):
    box = {}

    def _c_target(param):
        try:
            box["result"] = func(*args, **kwargs)
        except BaseException as e:
            box["error"] = e
        return 0

    c_callback = THREAD_ROUTINE(_c_target)
    thread_id = wintypes.DWORD()

    h_thread = kernel32.CreateThread(
        None,
        ctypes.c_size_t(STACK_SIZE_BYTES),
        c_callback,
        None,
        0,
        ctypes.byref(thread_id),
    )
    if not h_thread:
        raise ctypes.WinError(ctypes.get_last_error())

    INFINITE = 0xFFFFFFFF
    kernel32.WaitForSingleObject(h_thread, INFINITE)
    kernel32.CloseHandle(h_thread)

    if "error" in box:
        raise box["error"]
    if "result" not in box:
        raise RuntimeError("Worker thread died without returning a result "
                            "(still crashing even with a 1 GiB stack).")
    return box["result"]


with open(os.path.join(HERE, "grid_base_1m.pkl"), "rb") as f:
    base = pickle.load(f)
ny, nx = base["shape"]
res = base["res"]
STATUS = base["status"]
print(f"grid shape=({ny},{nx})  n_nodes={ny*nx}", flush=True)

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
CHANNEL_HEAD_M2 = 5.0e4


def _route_flow_worker(z_flat, tag):
    print(f"   [{tag}] building grid...", flush=True)
    grid = RasterModelGrid((ny, nx), xy_spacing=res)
    zf = grid.add_zeros("topographic__elevation", at="node")
    zf[:] = z_flat
    grid.status_at_node[:] = STATUS.copy()
    fa = FlowAccumulator(grid, flow_director="FlowDirectorD8",
                          depression_finder=DepressionFinderAndRouter)
    print(f"   [{tag}] routing flow (slow step)...", flush=True)
    fa.run_one_step()
    print(f"   [{tag}] done", flush=True)
    return grid


def route_flow(z_flat, tag=""):
    return run_in_big_stack(_route_flow_worker, z_flat, tag)


def fit_slope_area(da, slope, n_bins=18):
    mask = (da > CHANNEL_HEAD_M2) & (slope > 0) & np.isfinite(slope)
    if mask.sum() < 10:
        return np.nan, np.nan, None, None
    log_a = np.log10(da[mask])
    log_s = np.log10(slope[mask])
    bins = np.linspace(log_a.min(), log_a.max(), n_bins)
    idx = np.digitize(log_a, bins)
    bin_a, bin_s = [], []
    for i in range(1, len(bins)):
        sel = idx == i
        if sel.sum() >= 3:
            bin_a.append(np.median(log_a[sel]))
            bin_s.append(np.median(log_s[sel]))
    if len(bin_a) < 3:
        return np.nan, np.nan, log_a, log_s
    bin_a, bin_s = np.array(bin_a), np.array(bin_s)
    coeffs = np.polyfit(bin_a, bin_s, 1)
    return -coeffs[0], 10 ** coeffs[1], bin_a, bin_s


fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
results = {}
for model, ax in zip(["spim", "space"], axes):
    for k in KEYS_ORDERED:
        print(f"loading {model}_{k}.npz", flush=True)
        d = np.load(f"{RES_DIR}/{model}_{k}.npz")
        steps = sorted(int(s.split("_")[1]) for s in d.files if s.startswith("snap_"))
        zf = d[f"snap_{steps[-1]}"]
        grid = route_flow(zf, tag=f"{model}_{k}")
        da = grid.at_node["drainage_area"]
        slope = grid.at_node["topographic__steepest_slope"]
        theta_fit, ks_fit, bin_a, bin_s = fit_slope_area(da, slope)
        results[(model, k)] = (theta_fit, ks_fit)
        if bin_a is not None:
            ax.plot(bin_a, bin_s, "o-", color=COLORS[k], ms=4, lw=1.3,
                    label=rf"$K_{{br}}={K_VALUES[k]:.1e}$: $\theta$={theta_fit:.2f}")
    ax.set_xlabel(r"log$_{10}$ drainage area (m$^2$)")
    ax.set_ylabel(r"log$_{10}$ slope")
    ax.set_title(f"{model.upper()}: binned slope-area regression, t=1000 yr (1 m)", fontsize=10)
    ax.legend(fontsize=6.5, loc="lower left")

plt.tight_layout()
plt.savefig(f"{HERE}/fig_slope_area_1kyr.png", dpi=230, bbox_inches="tight")
print("Saved fig_slope_area_1kyr.png", flush=True)

print(f"\n=== Fitted theta (true = 0.5) and ks per K/model (1 m, 1000yr/DT=10yr) ===")
for (model, k), (theta_fit, ks_fit) in results.items():
    print(f"{model:6s} {k:14s} K_br={K_VALUES[k]:.2e}  theta_fit={theta_fit:.3f}  ks_fit={ks_fit:.3e}")
