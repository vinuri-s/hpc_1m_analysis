"""
1000yr/DT=10yr version of 07_channel_profiles.py: longitudinal profiles
and chi-plots for the trunk channel, at t=0 and t=1000 yr, across all 7 K
values and both models, for the results_1kyr_dt10/ production run.

Uses the same big-stack-thread workaround as the original 1 Myr script
(Windows native stack overflow on the full 1 m grid otherwise).

Usage:
    python3 07_channel_profiles_1kyr.py
"""
import os
import sys
import ctypes
from ctypes import wintypes
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from landlab import RasterModelGrid
from landlab.components import ChiFinder, FlowAccumulator, DepressionFinderAndRouter

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
OUTLET_ID = base["outlet_id"]
Z0 = base["z0"]
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
M_SP, N_SP = 0.5, 1.0
HEAD_THRESHOLD_M2 = 5.0e4


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
    print(f"   [{tag}] computing chi...", flush=True)
    cf = ChiFinder(grid, min_drainage_area=HEAD_THRESHOLD_M2, reference_concavity=M_SP / N_SP)
    cf.calculate_chi()
    print(f"   [{tag}] done", flush=True)
    return grid


def route_flow(z_flat, tag=""):
    return run_in_big_stack(_route_flow_worker, z_flat, tag)


def trace_trunk(grid):
    receiver = grid.at_node["flow__receiver_node"]
    da = grid.at_node["drainage_area"]
    n = grid.number_of_nodes
    donors_of = [[] for _ in range(n)]
    for node in range(n):
        r = receiver[node]
        if r != node:
            donors_of[r].append(node)
    path = [OUTLET_ID]
    current = OUTLET_ID
    while True:
        donors = donors_of[current]
        if not donors:
            break
        best = max(donors, key=lambda d: da[d])
        if da[best] < HEAD_THRESHOLD_M2:
            break
        path.append(best)
        current = best
    return np.array(path)


def profile_along_path(grid, path):
    xy = grid.xy_of_node
    x, y = xy[path, 0], xy[path, 1]
    dist = np.zeros(len(path))
    dist[1:] = np.cumsum(np.hypot(np.diff(x), np.diff(y)))
    z = grid.at_node["topographic__elevation"][path]
    chi = grid.at_node["channel__chi_index"][path]
    return dist, z, chi


fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex="col")

grid0 = route_flow(Z0, tag="t=0")
path0 = trace_trunk(grid0)
dist0, z0p, chi0 = profile_along_path(grid0, path0)
dist0 = dist0[-1] - dist0
order0 = np.argsort(dist0)
print("t=0 routing + trunk trace complete", flush=True)

for row, model in enumerate(["spim", "space"]):
    ax_prof, ax_chi = axes[row]
    ax_prof.plot(dist0[order0] / 1000, z0p[order0], color="k", ls="--", lw=1.3, label="t=0 (initial)")
    ax_chi.plot(chi0[order0], z0p[order0], color="k", ls="--", lw=1.3)
    for k in KEYS_ORDERED:
        print(f"loading {model}_{k}.npz", flush=True)
        d = np.load(f"{RES_DIR}/{model}_{k}.npz")
        steps = sorted(int(s.split("_")[1]) for s in d.files if s.startswith("snap_"))
        zf = d[f"snap_{steps[-1]}"]
        grid = route_flow(zf, tag=f"{model}_{k}")
        path = trace_trunk(grid)
        dist, z, chi = profile_along_path(grid, path)
        dist = dist[-1] - dist
        order = np.argsort(dist)
        ax_prof.plot(dist[order] / 1000, z[order], color=COLORS[k], lw=1.6,
                     label=rf"$K_{{br}}={K_VALUES[k]:.2e}$")
        ax_chi.plot(chi[order], z[order], color=COLORS[k], lw=1.6)
    ax_prof.set_ylabel(("SPIM\n" if row == 0 else "SPACE\n") + "Elevation (m)")
    ax_chi.set_ylabel("Elevation (m)")
    if row == 1:
        ax_prof.set_xlabel("Distance upstream from outlet (km)")
        ax_chi.set_xlabel(r"$\chi$ (m)")
    ax_prof.set_title(f"{model.upper()} longitudinal profile at t=1000 yr (1 m)", fontsize=10)
    ax_chi.set_title(rf"{model.upper()} $\chi$-elevation plot (1 m)", fontsize=10)

axes[0, 0].legend(fontsize=7, loc="lower right")
plt.tight_layout()
plt.savefig(f"{HERE}/fig_channel_profiles_1kyr.png", dpi=230, bbox_inches="tight")
print("Saved fig_channel_profiles_1kyr.png", flush=True)
