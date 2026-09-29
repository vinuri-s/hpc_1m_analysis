# 1 m native-resolution sensitivity sweep (HPC package)

This reproduces the paper's Section 4 sensitivity analysis (SPIM vs SPACE,
7 log-spaced basalt erodibility values from the Thor database) directly on
the **native 1 m LiDAR DEM**, with no downsampling. The version already in
the paper uses a 15 m block-mean-resampled grid (180x178 = 32,040 nodes)
because that's what fit inside a sandboxed session with per-call wall-clock
limits. This package is meant to run on your own machine/HPC to test
whether the paper's conclusions are resolution-dependent.

## What's identical to the paper's 15 m run

- K values (7, log-spaced): `3.024e-09, 1.191e-08, 4.693e-08, 1.849e-07, 7.283e-07, 2.869e-06, 1.130e-05` yr⁻¹, spanning UCS = 217.02 MPa (99th percentile) to 3.55 MPa (weakest measured) for basalt in the Thor database (n=341), converted via K ∝ UCS⁻², anchored to Quye-Sawyer et al. (2020). Note: the raw database maximum (976.81 MPa) is deliberately excluded -- it and the next-highest value (373.9 MPa) both trace to a single study (Teymen and Menguec 2020) at rebound numbers (R=99, R=84.6) at/near Schmidt-hammer saturation, isolated from the rest of the distribution and almost certainly a regression-extrapolation artifact rather than real rock strength. The 99th percentile sits right at the edge of the genuine, multi-source-corroborated cluster (~130-230 MPa, consistent with published basalt UCS ranges of 100-300 MPa).
- `m_sp = 0.5`, `n_sp = 1.0`
- `K_sed = 10 x K_br` in SPACE (per your original instruction)
- Uplift rate `U = 8.0e-6 m/yr` (calibrated from the real DEM's channel steepness at 15 m -- see note in `03_run_one_1m.py` if you'd rather recalibrate from the 1 m network)
- Single-outlet boundary conditions (all watershed-exterior + true grid-perimeter nodes closed except the lowest-elevation valid perimeter node)
- Total simulated time target: 1.0 Myr (the *timestep and number of steps used to reach it* are NOT fixed yet -- see "Required first step" below)

## What's different

- Grid: **2670 x 2700 = 7,209,000 total nodes** at 1 m, vs 32,040 at 15 m -- roughly **225x** more nodes. No block-mean resampling anywhere.
- Flow routing: **`FlowAccumulator` (plain D8) + `LakeMapperBarnes`** (fills pits and redirects/reaccumulates flow around them). This is NOT what the 15 m paper run uses (`DepressionFinderAndRouter`), and does NOT use richdem/`PriorityFloodFlowRouter`. Both of those were tried and rejected for this specific project:
  - **richdem**: on Windows, `pip install richdem` fails to compile (richdem 0.3.4 on PyPI bundles an old `pybind11` incompatible with modern Python's C API -- installing MSVC Build Tools does NOT fix this, it's a real source incompatibility). The conda-forge build installs and imports fine, but crashes at runtime (`IndexError: Out of bounds on buffer access`) specifically on this grid's boundary-condition layout.
  - **DepressionFinderAndRouter**: crashes natively at 1 m's node count (a raw Windows access-violation-style crash, not a catchable Python exception) -- it's not built for grids this large.
  - **LakeMapperBarnes** is part of landlab core (no extra dependency, nothing to install beyond `landlab` itself) and was verified to run cleanly through 100 test steps on the full 7.2M-node grid with no crash and no NaNs.

## Required first step: benchmark before committing to a real run

**Do not just edit `DT`/`N_STEPS_TARGET` in `03_run_one_1m.py` and launch the full sweep.** Raw flow-routing at 7.2M nodes is slow -- in testing, a single `FlowAccumulator` step alone (before any depression handling) took over 30 seconds. At the original plan (2000 steps of 500 yr), that's 16+ hours *per config*, and with 14 configs run sequentially on one machine, that's over a week.

Run `test_lakemapper_1m.py` first:
```
python test_lakemapper_1m.py
```
This runs 100 short (10 yr) steps on the full 1 m grid using the same `FlowAccumulator`+`LakeMapperBarnes` routing as the real runner, prints per-step timing as it goes, and at the end prints extrapolated total-runtime estimates for both a larger-timestep plan (`DT=5000, N_STEPS_TARGET=200`, same 1 Myr total, 10x fewer routing calls) and the original plan (`DT=500, N_STEPS_TARGET=2000`). **Use those real numbers, not any estimate written in this README or in code comments, to decide what `DT`/`N_STEPS_TARGET` to actually use** -- edit the placeholder values in `03_run_one_1m.py` (clearly marked) once you know.

Note: `StreamPowerEroder` and `SpaceLargeScaleEroder` both use implicit/semi-implicit solvers that remain numerically stable at large timesteps, so increasing `DT` well above 500 yr is a legitimate way to cut runtime, not a hack -- the trade-off is coarser temporal resolution (less precision on exactly when transitions between response regimes occur), not instability.

## Computational-cost warning

- **Runtime**: depends entirely on what `test_lakemapper_1m.py` measures on your hardware -- see above. Don't guess.
- **Memory**: each `float64` node field is ~58 MB at this grid size. Landlab keeps several such fields per component (elevation, slope, receiver, drainage area, and for SPACE also soil depth, bedrock elevation, sediment flux, etc.) -- budget several GB per running config, more for SPACE than SPIM.
- **Disk**: each snapshot saved is another ~58 MB; with `SNAPSHOT_EVERY=200` (or whatever you set it to after the benchmark) that's several snapshots/config, several hundred MB/config for SPIM and roughly double for SPACE (also saves soil depth) -- full output across all 14 configs could reach several tens of GB. Increase `SNAPSHOT_EVERY` in `03_run_one_1m.py` before running if disk is tight.
- **OneDrive**: if this folder lives inside a OneDrive-synced path, pause OneDrive sync before running (Settings tray icon -> Pause syncing). Live syncing while scripts write large binary files (checkpoints, `.npz` results) has caused intermittent "file not found" errors in testing.

## How to run

1. `watershed_of_interest.tif` and `Search tool.xlsx` are already included in this folder.
2. Create the conda environment (no richdem needed):
   ```
   conda create -n erodibility python=3.11 -y
   conda activate erodibility
   conda install -c conda-forge landlab numpy scipy rasterio matplotlib pandas openpyxl -y
   ```
   If `conda activate erodibility` errors with something about PowerShell not being recognized, run `conda init cmd.exe` once, then close and reopen Anaconda Prompt and try again. If you instead see `(erodibility)` appear in your prompt alongside a PowerShell error, that error is harmless noise (activation still worked) -- confirm with `where python` (should show a path containing `envs\erodibility`).
3. Run setup once (cheap): `python 01_setup_grid_1m.py` then `python 02_build_grid_1m.py`. Check the printed diagnostics: exactly one `FIXED_VALUE` outlet node, drainage area close to the full watershed extent, and `Flow router used: FlowAccumulator+LakeMapperBarnes` with no crash/traceback.
4. Run `test_lakemapper_1m.py` and use its output to set `DT`/`N_STEPS_TARGET` in `03_run_one_1m.py` (see "Required first step" above). Don't skip this.
5. Run the full sweep with `run_all.py` (plain Python, works the same on Windows/Mac/Linux -- prefer this over `run_all_windows.bat`, which has a known bug in its retry-loop logic and shouldn't be used):
   ```
   python run_all.py
   ```
   It skips setup steps that already succeeded, runs all 14 model/K configs (resuming from checkpoint if interrupted), then all the figure scripts, and stops with a clear message if anything actually fails rather than looping uselessly.
   - If you have several free CPU cores and want to speed things up by running configs in parallel instead of sequentially, open multiple Anaconda Prompt windows (each with `conda activate erodibility`) and run individual `python 03_run_one_1m.py <model> <kname>` calls directly, one per window, instead of `run_all.py`.
   - If you have access to a Linux HPC with SLURM (in addition to or instead of this Windows machine), `submit_array.slurm` is available, but will need the same richdem-removal edits applied to work (it currently assumes the old routing setup) -- ask if you want this updated too.
6. Once all 14 `results/{model}_{kname}.npz` + `_log.json` files exist (already done if you used `run_all.py`), the figure scripts run in this order:
   - `python 04_make_figures_1m.py` -- core comparison figures (timeseries, diffmaps by K, structural diffmaps)
   - `python 05_ucs_distribution.py` -- Thor basalt UCS histogram + K sampling/outlier-exclusion figure (doesn't depend on the sweep, can be run any time)
   - `python 06_sediment_diagnostics.py` -- SPACE soil-depth maps, sediment-flux/soil-depth time series, mass-balance sanity check
   - `python 07_channel_profiles.py` -- longitudinal profile + chi-plots for the trunk channel (re-routes flow on each snapshot; the trunk-tracing step uses a per-node Python loop, noticeably slower at 1 m than at 15 m)
   - `python 08_final_state_maps.py` -- hillshaded final-elevation maps (actual topography, not differences)
   - `python 09_slope_area.py` -- binned slope-area regression with fitted concavity/steepness index per K

## What to send back

The `results/` folder (`.npz` and `_log.json` files -- `ckpt_*.pkl` checkpoints aren't needed once a config is complete and can be deleted) plus all generated `fig_*_1m.png`/`.pdf` files. Also mention what `DT`/`N_STEPS_TARGET` you ended up using (from the benchmark step) so the comparison against the 15 m results accounts for any difference in temporal resolution. I'll use these to check whether the 15 m paper results (three-regime response: K1-K4 insensitive, K5 transitional, K6 strong divergence with non-monotonic SPACE relief, K7 near-total base-levelling) hold up at native 1 m resolution, or whether the 15 m grid was smoothing out something resolution-dependent.
