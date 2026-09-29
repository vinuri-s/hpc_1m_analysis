"""
Checkpointable single-configuration sensitivity runner -- native 1 m
resolution, SHORT-DURATION variant.

Same physics, K values, m/n exponents, uplift rate, and grid as
03_run_one_1m.py (the main 1 Myr / DT=5000yr production run). The only
difference is the timestep and total duration:

    DT = 10 yr, N_STEPS_TARGET = 100  ->  1,000 yr total model time

This is a MUCH finer timestep (500x smaller than the production run's
DT=5000) and a MUCH shorter total duration (1000x shorter than the
production run's 1 Myr). Two things to keep in mind about this combo:

  1. Finer DT should make SPACE more numerically stable (less likely to
     show the K6/K7 pit-filling blowup seen in the production run) --
     this is the main thing worth checking with this variant.
  2. At only 1000 yr with uplift_rate=8.0e-6 m/yr, total imposed uplift
     is just 0.008 m over the whole run, and this catchment's K-dependent
     divergence in the production run only became visually/numerically
     obvious after several hundred thousand years. So relief/slope/
     mean-erosion-rate differences between K values will likely be small
     to negligible at this duration -- useful as a quick stability check,
     not a substitute for the full 1 Myr sweep.

Writes to its own results_1kyr_dt10/ subfolder -- does NOT touch
results/ (1 Myr/DT=5000) or results_10kyr_dt10/ (10,000yr/DT=10) from
earlier variants.

Usage:
    python3 03_run_one_1m_1kyr.py <model: spim|space> <kname> [wall_budget_s]
"""
import sys
import os
import json
import time
import pickle
import threading
import traceback
import numpy as np

from landlab import RasterModelGrid
from landlab.components import (
    FlowAccumulator,
    LakeMapperBarnes,
    StreamPowerEroder,
    SpaceLargeScaleEroder,
    ExponentialWeatherer,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results_1kyr_dt10")
os.makedirs(RESULTS_DIR, exist_ok=True)

with open(os.path.join(HERE, "grid_base_1m.pkl"), "rb") as f:
    BASE = pickle.load(f)

ny, nx = BASE["shape"]
RES = BASE["res"]
VALID = BASE["valid_flat"]
STATUS = BASE["status"]
Z0 = BASE["z0"]
OUTLET_ID = BASE["outlet_id"]

# Identical to the main 1m sweep: 7 log-spaced K values spanning the
# Thor-database (Haag & Schoenbohm 2025) basalt UCS range (n=341),
# converted via K = K_ref*(UCS_ref/UCS)^2 anchored at K_ref=1.0e-7 yr^-1,
# UCS_ref=37.74 MPa (Quye-Sawyer et al. 2020 Sardinia basalt).
K_VALUES = {
    "K1_3.02e-09": 3.024e-09,   # UCS = 217.0 MPa (99th-pct basalt UCS)
    "K2_1.19e-08": 1.191e-08,
    "K3_4.69e-08": 4.693e-08,
    "K4_1.85e-07": 1.849e-07,
    "K5_7.28e-07": 7.283e-07,
    "K6_2.87e-06": 2.869e-06,
    "K7_1.13e-05": 1.130e-05,   # UCS = 3.55 MPa (weakest basalt measured)
}

M_SP = 0.5
N_SP = 1.0
UPLIFT_RATE = 8.0e-6  # m/yr -- same physical rate as all other variants

# --- The only real change vs. 03_run_one_1m.py ---
DT = 10.0              # yr  (500x finer than the 1 Myr production run)
N_STEPS_TARGET = 100    # -> 1,000 yr total model time
SNAPSHOT_EVERY = 10     # -> 11 snapshots across the run
LOG_EVERY = 1
CKPT_EVERY = 10

# At ~36 s/step (routing-dominated, same cost regardless of DT), 100
# steps is ~1 hour per config, ~14 hours for all 14 sequentially, or a
# few hours if you run several in parallel (see
# run_all_1kyr_parallel.ps1).


def checkpoint_path(model, kname):
    return os.path.join(RESULTS_DIR, f"ckpt_{model}_{kname}.pkl")


def make_flow_router(grid):
    """FlowAccumulator (plain D8) + LakeMapperBarnes -- the only router
    validated to run cleanly on the full 7.2M-node 1 m grid on this
    Windows/landlab build (PriorityFloodFlowRouter fails to build/crashes;
    DepressionFinderAndRouter causes a native Windows stack overflow at
    this node count)."""
    fa = FlowAccumulator(grid, flow_director="FlowDirectorD8")
    lmb = LakeMapperBarnes(grid, method="Steepest", fill_flat=False,
                            redirect_flow_steepest_descent=True,
                            reaccumulate_flow=True, track_lakes=True,
                            ignore_overfill=True)
    return fa, lmb, "FlowAccumulator+LakeMapperBarnes"


def route_flow(fa, lmb):
    fa.run_one_step()
    lmb.run_one_step()


def make_components(grid, model, K):
    fa, lmb, router_name = make_flow_router(grid)
    if model == "spim":
        sp = StreamPowerEroder(grid, K_sp=K, m_sp=M_SP, n_sp=N_SP)
        return fa, lmb, sp, None, router_name
    else:
        ew = ExponentialWeatherer(grid, soil_production_maximum_rate=2.0e-4,
                                   soil_production_decay_depth=0.5)
        sp = SpaceLargeScaleEroder(grid, K_sed=10.0 * K, K_br=K, F_f=0.0, phi=0.0,
                                    H_star=1.0, v_s=1.0, m_sp=M_SP, n_sp=N_SP,
                                    sp_crit_sed=0.0, sp_crit_br=0.0)
        return fa, lmb, sp, ew, router_name


def log_point(state, grid, z, core, model):
    zc = z[VALID]
    slope = grid.at_node["topographic__steepest_slope"][VALID]
    dt_elapsed = DT * LOG_EVERY if state["step"] > 0 else DT
    erosion_rate = float(np.mean((state["z_prev"][VALID] - z[VALID]) / dt_elapsed))
    log = state["log"]
    log["t"].append(state["step"] * DT)
    log["relief"].append(float(np.nanmax(zc) - np.nanmin(zc)))
    log["mean_slope"].append(float(np.nanmean(slope)))
    log["median_slope"].append(float(np.nanmedian(slope)))
    log["mean_elev_core"].append(float(np.nanmean(z[core])))
    log["mean_erosion_rate"].append(erosion_rate)
    if model == "space":
        soil = grid.at_node["soil__depth"]
        qs_field = grid.at_node["sediment__influx"]
        log["mean_soil_depth"].append(float(np.nanmean(soil[VALID])))
        log["sediment_flux_outlet"].append(float(qs_field[OUTLET_ID]))
    state["z_prev"] = z.copy()


def save_checkpoint(ckpt_file, z, soil, bedrock, state, model):
    ck_out = {"z": np.array(z), "state": state}
    if model == "space":
        ck_out["soil"] = np.array(soil)
        ck_out["bedrock"] = np.array(bedrock)
    with open(ckpt_file, "wb") as f:
        pickle.dump(ck_out, f)


def main():
    model = sys.argv[1]
    kname = sys.argv[2]
    wall_budget = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0e12
    K = K_VALUES[kname]

    t_start = time.time()
    ckpt_file = checkpoint_path(model, kname)

    grid = RasterModelGrid((ny, nx), xy_spacing=RES)
    z = grid.add_zeros("topographic__elevation", at="node")
    grid.status_at_node[:] = STATUS.copy()
    core = grid.status_at_node == grid.BC_NODE_IS_CORE

    soil = bedrock = None
    if model == "space":
        soil = grid.add_zeros("soil__depth", at="node")
        bedrock = grid.add_zeros("bedrock__elevation", at="node")

    if os.path.exists(ckpt_file):
        with open(ckpt_file, "rb") as f:
            ck = pickle.load(f)
        z[:] = ck["z"]
        if model == "space":
            soil[:] = ck["soil"]
            bedrock[:] = ck["bedrock"]
        state = ck["state"]
        print(f"Resumed {model}/{kname} from step {state['step']}", flush=True)
    else:
        z[:] = Z0.copy()
        if model == "space":
            soil[VALID] = 1.0
            bedrock[:] = z - soil
        state = {"step": 0, "log": {k: [] for k in
                  ["t", "relief", "mean_slope", "median_slope", "mean_elev_core",
                   "mean_erosion_rate"] + (["mean_soil_depth", "sediment_flux_outlet"] if model == "space" else [])},
                  "snapshots": {}, "snapshots_soil": {} if model == "space" else None,
                  "z_prev": z.copy()}
        print(f"Starting fresh {model}/{kname} at 1 m native resolution, "
              f"1000yr/DT=10yr variant ({ny}x{nx} = {ny*nx:,} nodes)", flush=True)

    fa, lmb, sp, ew, router_name = make_components(grid, model, K)
    print("Using flow router:", router_name, flush=True)
    route_flow(fa, lmb)

    if state["step"] == 0:
        log_point(state, grid, z, core, model)
        state["snapshots"][0] = z.copy()
        if model == "space":
            state["snapshots_soil"][0] = soil.copy()

    last_ckpt_step = state["step"]
    while state["step"] < N_STEPS_TARGET:
        if time.time() - t_start > wall_budget:
            break
        state["step"] += 1
        if model == "spim":
            z[core] += UPLIFT_RATE * DT
            route_flow(fa, lmb)
            sp.run_one_step(dt=DT)
        else:
            bedrock[core] += UPLIFT_RATE * DT
            z[:] = bedrock + soil
            route_flow(fa, lmb)
            ew.calc_soil_prod_rate()
            prod = grid.at_node["soil_production__rate"]
            soil[core] += prod[core] * DT
            bedrock[core] -= prod[core] * DT
            sp.run_one_step(dt=DT)

        if state["step"] % LOG_EVERY == 0 or state["step"] == N_STEPS_TARGET:
            log_point(state, grid, z, core, model)

        if state["step"] % SNAPSHOT_EVERY == 0 or state["step"] == N_STEPS_TARGET:
            state["snapshots"][state["step"]] = z.copy()
            if model == "space":
                state["snapshots_soil"][state["step"]] = soil.copy()

        if state["step"] - last_ckpt_step >= CKPT_EVERY:
            save_checkpoint(ckpt_file, z, soil, bedrock, state, model)
            last_ckpt_step = state["step"]
            print(f"  [safety checkpoint] step {state['step']}/{N_STEPS_TARGET}", flush=True)

    save_checkpoint(ckpt_file, z, soil, bedrock, state, model)

    elapsed = time.time() - t_start
    done = state["step"] >= N_STEPS_TARGET
    print(f"{model}/{kname}: step {state['step']}/{N_STEPS_TARGET} "
          f"({elapsed:.1f}s this call)", flush=True)

    if done:
        snaps = state["snapshots"]
        savez_kwargs = {f"snap_{k}": v for k, v in snaps.items()}
        if model == "space":
            savez_kwargs.update({f"soil_{k}": v for k, v in state["snapshots_soil"].items()})
        np.savez(os.path.join(RESULTS_DIR, f"{model}_{kname}.npz"), **savez_kwargs)
        with open(os.path.join(RESULTS_DIR, f"{model}_{kname}_log.json"), "w") as f:
            json.dump(state["log"], f)
        print("STATUS: COMPLETE", flush=True)
    else:
        print("STATUS: INCOMPLETE (call again with the same arguments to resume)", flush=True)


# --- Windows native-crash workaround (same as 03_run_one_1m.py) --------
STACK_SIZE_BYTES = 1024 * 1024 * 1024  # 1 GiB -- ctypes CreateThread, not
                                        # threading.stack_size(), since that
                                        # API rejects all nonzero sizes on
                                        # this Windows/Python build.

import ctypes
from ctypes import wintypes

THREAD_ROUTINE = ctypes.WINFUNCTYPE(wintypes.DWORD, wintypes.LPVOID)


def _run_in_big_stack_thread():
    error_box = {}

    def _target(_):
        try:
            main()
        except Exception:
            error_box["exc"] = traceback.format_exc()
        return 0

    callback = THREAD_ROUTINE(_target)
    thread_id = wintypes.DWORD()
    handle = ctypes.windll.kernel32.CreateThread(
        None, ctypes.c_size_t(STACK_SIZE_BYTES), callback, None, 0,
        ctypes.byref(thread_id)
    )
    if not handle:
        raise ctypes.WinError(ctypes.get_last_error())

    ctypes.windll.kernel32.WaitForSingleObject(handle, 0xFFFFFFFF)
    ctypes.windll.kernel32.CloseHandle(handle)

    if "exc" in error_box:
        print("\n=== Exception inside worker thread ===", file=sys.stderr)
        print(error_box["exc"], file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _run_in_big_stack_thread()
