# Native 1 m resolution sensitivity analysis

Code for the native-resolution (1 m) landscape evolution sensitivity
analysis reported in Section 4 of:

> Piyathilake, V., Hughes, M. W., and Wilson, M.: Quantifying Erodibility
> in Stream Power Incision Models: Practical Guidance and Model Response
> Dynamics, *Earth Surface Dynamics* (in review).

The analysis compares two landscape evolution model formulations —
detachment-limited (SPIM, `StreamPowerEroder`) and sediment-flux-dependent
(SPACE, `SpaceLargeScaleEroder`), both from
[Landlab](https://landlab.readthedocs.io/) — run on a real 1 m LiDAR
digital elevation model (DEM) of a basalt catchment near Little River,
Banks Peninsula, New Zealand, across seven log-spaced bedrock erodibility
($K_\mathrm{br}$) values spanning the strength range of basalt in the
Thor rock-strength database (Haag et al., 2025). Fourteen
model/$K_\mathrm{br}$ configurations are each run for 1000 years of
model time at a 10-year timestep, at the DEM's native 1 m grid resolution
(no spatial resampling), directly reproducing the figures and summary
tables in Section 4 and its supplement.

## What this reproduces

| Script | Output figure(s) | Figure number | Location |
|---|---|---|---|
| `05_ucs_distribution.py` | `fig_ucs_distribution_1m.png/.pdf` | Figure 4 | Main paper |
| `04_make_figures_1kyr.py` | `fig_sensitivity_timeseries_1kyr.png` | Figure 5 | Main paper |
| `04_make_figures_1kyr.py` | `fig_diffmaps_by_K_1kyr.png` | Figure 6 | Main paper |
| `04_make_figures_1kyr.py` | `fig_diffmaps_structural_1kyr.png` | Figure 7 | Main paper |
| `06_sediment_diagnostics_1kyr.py` | `fig_soil_depth_maps_1kyr.png` | Figure S6 | Supplement |
| `06_sediment_diagnostics_1kyr.py` | `fig_sediment_flux_soil_ts_1kyr.png` | Figure S7 | Supplement |
| `07_channel_profiles_1kyr.py` | `fig_channel_profiles_1kyr.png` | Figure S8 | Supplement |
| `09_slope_area_1kyr.py` | `fig_slope_area_1kyr.png` | Figure S9 | Supplement |
| `08_final_state_maps_1kyr.py` | `fig_final_state_maps_1kyr.png` | Figure S10 | Supplement |

(Figure and table numbers follow the current manuscript draft and may
shift if content is added, removed, or reordered before publication.)

Three tables are also compiled from the `results_1kyr_dt10/*_log.json`
files written by `03_run_one_1m_1kyr.py`:

| Content | Table number | Location |
|---|---|---|
| Summary of all 14 sensitivity configurations after 1000 yr (relief, mean slope, mean erosion rate, eroded volume) | Table 4 | Main paper |
| Total SPACE eroded volume per $K_\mathrm{br}$ | Table 5 | Main paper |
| SPACE mass-balance check (eroded volume vs. cumulative sediment export) | Table S2 | Supplement |

Other figures in the paper (catchment-overview and location maps, the
conceptual/schematic diagrams, and the coarser 15 m sensitivity run used
elsewhere in the manuscript) are produced outside this package and are
not covered here.

## Model configuration

- Seven log-spaced $K_\mathrm{br}$ values (yr$^{-1}$): `3.024e-09,
  1.191e-08, 4.693e-08, 1.849e-07, 7.283e-07, 2.869e-06, 1.130e-05`,
  spanning basalt UCS = 217.0 MPa (99th percentile of n=341 Thor
  measurements, used as the resistant anchor) to 3.55 MPa (database
  minimum, erodible anchor), converted via $K \propto \mathrm{UCS}^{-2}$
  and anchored to the independent estimate of Quye-Sawyer et al. (2020).
  Two Schmidt-hammer-saturation outliers (976.8 and 373.9 MPa, both from
  Teymen and Mengüç, 2020) are excluded as regression-extrapolation
  artifacts rather than genuine rock-strength measurements; see
  `05_ucs_distribution.py` for the full derivation.
- SPACE sediment erodibility $K_\mathrm{sed} = 10 \times K_\mathrm{br}$.
- $m_\mathrm{sp} = 0.5$, $n_\mathrm{sp} = 1.0$.
- Uplift rate $U = 8.0\times10^{-6}$ m yr$^{-1}$.
- Single-outlet boundary conditions: all watershed-exterior nodes and the
  true grid perimeter are closed except the lowest-elevation valid
  perimeter node (fixed-value outlet).
- Flow routing: `FlowAccumulator` (D8) + `LakeMapperBarnes` for
  depression filling/redirection. `LakeMapperBarnes` was selected over
  `DepressionFinderAndRouter` and richdem's `PriorityFloodFlowRouter`
  because both failed to run reliably at this grid's node count
  (~7.2 million nodes); `LakeMapperBarnes` is Landlab-native and was
  verified stable on the full grid.
- Timestep: `DT = 10` yr, `N_STEPS_TARGET = 100` (1000 yr total model
  time).

## Repository layout

```
01_setup_grid_1m.py            Load the 1 m DEM, build the valid-cell mask
02_build_grid_1m.py             Build the Landlab grid, boundary conditions,
                                 initial flow routing
03_run_one_1m_1kyr.py           Run one (model, K) configuration; checkpointable
04_make_figures_1kyr.py         Timeseries + elevation-difference figures
05_ucs_distribution.py          Thor UCS distribution / K-sampling figure
06_sediment_diagnostics_1kyr.py Soil-depth maps, sediment-flux time series
07_channel_profiles_1kyr.py     Longitudinal profile + chi-plots
08_final_state_maps_1kyr.py     Hillshaded final-elevation maps
09_slope_area_1kyr.py           Binned slope-area regression per K
run_all_1kyr_auto.ps1           Windows launcher for all 14 configs (memory-throttled)
requirements.txt                Python dependencies
```

## Input data

- **`watershed_of_interest.tif`** (included in this repository) — the 1 m
  LiDAR DEM of the case-study catchment, clipped to the watershed extent.
  See "Data licensing and attribution" below for provenance and licence
  terms.
- **`Search tool.xlsx`** (not included — place in the repository root
  before running `05_ucs_distribution.py`) — the basalt subset of the
  Thor rock-strength database (Haag et al., 2025), providing the 341 UCS
  measurements used to derive the sampled $K_\mathrm{br}$ range.

## Requirements

Python 3.11, with dependencies listed in `requirements.txt`:

```bash
conda create -n erodibility python=3.11 -y
conda activate erodibility
conda install -c conda-forge landlab numpy scipy rasterio matplotlib pandas openpyxl -y
```

richdem is deliberately not used (see "Model configuration" above and
comments in `02_build_grid_1m.py`).

## Running the pipeline

Place `Search tool.xlsx` in the repository root (`watershed_of_interest.tif`
is already included), then run:

```bash
python 01_setup_grid_1m.py
python 02_build_grid_1m.py
```

Run all fourteen (model, $K_\mathrm{br}$) configurations. Each takes the
model (`spim` or `space`) and a K label as arguments and is
checkpointable/resumable:

```bash
python 03_run_one_1m_1kyr.py spim K1_3.02e-09
python 03_run_one_1m_1kyr.py spim K2_1.19e-08
# ... through K7_1.13e-05, then the same seven for "space"
```

On Windows, `run_all_1kyr_auto.ps1` runs all fourteen configurations
automatically, throttled to stay within available system memory:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\run_all_1kyr_auto.ps1
```

Once all fourteen `results_1kyr_dt10/{model}_{kname}.npz` and
`_log.json` files exist, generate the figures (order does not matter,
except that `05_ucs_distribution.py` has no dependency on the model runs
and can be run at any time):

```bash
python 04_make_figures_1kyr.py
python 05_ucs_distribution.py
python 06_sediment_diagnostics_1kyr.py
python 07_channel_profiles_1kyr.py
python 08_final_state_maps_1kyr.py
python 09_slope_area_1kyr.py
```

## Computational cost

- **Grid size**: 2670 x 2700 = 7,209,000 total nodes at 1 m resolution.
- **Memory**: each `float64` node field is ~58 MB at this grid size;
  Landlab keeps several fields per running configuration (elevation,
  slope, receiver, drainage area, and for SPACE also soil depth, bedrock
  elevation, and sediment flux) — budget several GB per configuration,
  more for SPACE than SPIM.
- **Disk**: each saved snapshot is ~58 MB; full output across all
  fourteen configurations can reach the tens-of-GB range depending on
  the configured snapshot interval.

## Citation

If you use this code, please cite the paper above. Please also cite:

- **Landlab**, the modelling toolkit this code is built on:
  - Hobley, D. E. J., Adams, J. M., Nudurupati, S. S., Hutton, E. W. H.,
    Gasparini, N. M., Istanbulluoglu, E., and Tucker, G. E.: Creative
    computing with Landlab: an open-source toolkit for building,
    coupling, and exploring two-dimensional numerical models of
    Earth-surface dynamics, *Earth Surface Dynamics*, 5, 21–46,
    <https://doi.org/10.5194/esurf-5-21-2017>, 2017.
  - Barnhart, K. R., Hutton, E. W. H., Tucker, G. E., Gasparini, N. M.,
    Istanbulluoglu, E., Hobley, D. E. J., Lyons, N. J., Mouchene, M.,
    Nudurupati, S. S., Adams, J. M., and Bandaragoda, C.: Short
    communication: Landlab v2.0: a software package for Earth surface
    dynamics, *Earth Surface Dynamics*, 8, 379–397,
    <https://doi.org/10.5194/esurf-8-379-2020>, 2020.
- **SPACE**, the sediment-flux-dependent erosion component used here:
  - Shobe, C. M., Tucker, G. E., and Barnhart, K. R.: The SPACE 1.0
    model: a Landlab component for 2-D calculation of sediment
    transport, bedrock erosion, and landscape evolution, *Geoscientific
    Model Development*, 10, 4577–4604,
    <https://doi.org/10.5194/gmd-10-4577-2017>, 2017.
- The underlying data sources: the LINZ 1 m LiDAR DEM, the NIWA REC2
  watersheds layer, and the Thor rock-strength database (Haag et al.,
  2025) — see "Data licensing and attribution" below for the DEM and
  watershed-boundary citations.

## License

The **code** in this repository (all `.py` and `.ps1` files) is licensed
under the MIT License — see `LICENSE`.

The **DEM** (`watershed_of_interest.tif`) is third-party data and is
**not** covered by the MIT licence above. It remains licensed under its
original terms, CC BY 4.0, as described below. If you reuse or
redistribute it, attribution to the original data providers is required
regardless of how the code is licensed.

## Data licensing and attribution

**`watershed_of_interest.tif`** (1 m LiDAR DEM, clipped to the
case-study catchment):

> The file `watershed_of_interest.tif` is derived from the New Zealand
> LiDAR 1 m DEM (tile BY24), sourced from the LINZ Data Service
> (<https://data.linz.govt.nz/layer/121859-new-zealand-lidar-1m-dem/>)
> and licensed for reuse under CC BY 4.0
> (<https://creativecommons.org/licenses/by/4.0/>). It has been clipped
> to the case-study catchment boundary. This redistribution is not
> endorsed by LINZ.

The catchment boundary used to clip the DEM is itself third-party data:

> The catchment boundary (watershed_1) is derived from the New Zealand
> Watersheds layer of NIWA's River Environment Classification version 2
> (REC2 v2.5; NIWA, Taihoro Nukurangi, 2018), available at
> <https://data-niwa.opendata.arcgis.com/datasets/watersheds> and
> licensed for reuse under CC BY 4.0
> (<https://creativecommons.org/licenses/by/4.0/>). The catchment of
> interest was extracted from this layer. This redistribution is not
> endorsed by NIWA.

**`Search tool.xlsx`** (not included — see "Input data" above): the
basalt subset of the Thor rock-strength database. Cite Haag et al.
(2025) if you obtain and use this file.
