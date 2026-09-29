"""
Step 2 (native-resolution version): build the Landlab RasterModelGrid
directly from the 1 m DEM (no resampling), set boundary conditions (closed
outside the watershed mask, single open outlet at the lowest-elevation
valid perimeter cell), run one flow-routing pass, and save the grid +
flow fields for the sensitivity runner.

Identical boundary-condition logic to the 15 m version (sens_02_grid.py),
just at native resolution. Uses FlowAccumulator (basic D8, no depression
handling) followed by LakeMapperBarnes to fill/redirect around pits.
richdem's PriorityFloodFlowRouter was tried and rejected for this
project: on Windows it either fails to build from source (old pybind11 in
richdem 0.3.4 incompatible with modern Python) or, when installed via
conda, crashes at runtime on this grid's boundary-condition layout
(IndexError in the underlying C++ code). DepressionFinderAndRouter (the
other pure-Python option) was also rejected: it crashes natively
(Windows-level access violation, not a catchable Python exception) at
this node count. LakeMapperBarnes is landlab-native (no external
dependency) and was verified to run cleanly on the full 1 m grid.

Usage:
    python3 02_build_grid_1m.py
"""
import os
import numpy as np
import pickle

from landlab import RasterModelGrid

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = HERE

data = np.load(os.path.join(OUT_DIR, "dem_1m.npz"))
z = data["z"]
valid = data["valid"]
res = float(data["res"])
ny, nx = z.shape
print("Grid shape:", ny, nx, "res:", res, "-> total nodes:", ny * nx)

grid = RasterModelGrid((ny, nx), xy_spacing=res)
zf = grid.add_zeros("topographic__elevation", at="node")

z_flat = z.flatten()
valid_flat = valid.flatten()
fill_value = np.nanmin(z) - 1.0
zf[:] = np.where(np.isfinite(z_flat), z_flat, fill_value)

# --- Boundary conditions (identical logic to the 15 m grid) --------------
grid.status_at_node[:] = grid.BC_NODE_IS_CORE
grid.status_at_node[~valid_flat] = grid.BC_NODE_IS_CLOSED
grid.status_at_node[grid.boundary_nodes] = grid.BC_NODE_IS_CLOSED

valid_ids = np.where(valid_flat)[0]
neighbours = grid.adjacent_nodes_at_node
is_edge = np.zeros(grid.number_of_nodes, dtype=bool)
for nid in valid_ids:
    for nb in neighbours[nid]:
        if nb == -1 or not valid_flat[nb]:
            is_edge[nid] = True
            break

edge_valid_ids = valid_ids[is_edge[valid_ids]]
outlet_id = edge_valid_ids[np.argmin(zf[edge_valid_ids])]
print("Outlet node id:", outlet_id, "elevation:", zf[outlet_id],
      "row,col:", outlet_id // nx, outlet_id % nx)

grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE

n_open = np.sum(grid.status_at_node == grid.BC_NODE_IS_FIXED_VALUE)
n_closed = np.sum(grid.status_at_node == grid.BC_NODE_IS_CLOSED)
n_core = np.sum(grid.status_at_node == grid.BC_NODE_IS_CORE)
print("n_open(fixed-value):", n_open, "n_closed:", n_closed, "n_core:", n_core)
assert n_open == 1, "Expected exactly one open (FIXED_VALUE) outlet node."

# --- Flow routing: FlowAccumulator (plain D8, no depression handling)
#     followed by LakeMapperBarnes to fill pits and redirect/reaccumulate
#     flow around them. See module docstring for why this replaced
#     richdem/DepressionFinderAndRouter. -----------------------------------
from landlab.components import FlowAccumulator, LakeMapperBarnes

fa = FlowAccumulator(grid, flow_director="FlowDirectorD8")
fa.run_one_step()
lmb = LakeMapperBarnes(grid, method="Steepest", fill_flat=False,
                        redirect_flow_steepest_descent=True,
                        reaccumulate_flow=True, track_lakes=True,
                        ignore_overfill=True)
lmb.run_one_step()
router_used = "FlowAccumulator + LakeMapperBarnes"

print("Flow router used:", router_used)

da = grid.at_node["drainage_area"]
print("Max drainage area (m^2):", np.nanmax(da), "-> km^2:", np.nanmax(da) / 1e6)
print("Total valid cell area (km^2):", valid.sum() * res * res / 1e6)

slope = grid.at_node["topographic__steepest_slope"]
print("Mean slope (valid nodes):", np.nanmean(slope[valid_flat]))

with open(os.path.join(OUT_DIR, "grid_base_1m.pkl"), "wb") as f:
    pickle.dump(
        {
            "z0": np.array(zf),
            "status": np.array(grid.status_at_node),
            "shape": (ny, nx),
            "res": res,
            "outlet_id": int(outlet_id),
            "valid_flat": np.array(valid_flat),
            "drainage_area": np.array(da),
            "slope": np.array(slope),
            "router_used": router_used,
        },
        f,
    )
print("Saved grid_base_1m.pkl")
