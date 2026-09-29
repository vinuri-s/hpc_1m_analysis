"""
Step 1 (native-resolution version): load the 1 m LiDAR DEM AS-IS -- no
block-mean resampling -- and save an .npz with elevation + valid-watershed
mask for the grid-building step.

This is the no-downsampling counterpart of sens_01_setup.py used for the
15 m sensitivity sweep in the paper. Everything else (K values, m, n,
uplift rate, timestep, number of steps) is unchanged, so results from this
1 m run are directly comparable to the 15 m results already in the paper.

Usage:
    python3 01_setup_grid_1m.py

Expects watershed_of_interest.tif in the same directory as this script
(or edit DEM_PATH below).
"""
import os
import numpy as np
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
DEM_PATH = os.path.join(HERE, "watershed_of_interest.tif")
OUT_DIR = HERE
os.makedirs(OUT_DIR, exist_ok=True)

with rasterio.open(DEM_PATH) as src:
    z = src.read(1)
    nodata = src.nodata
    res = src.res[0]
    print("DEM shape:", z.shape, "res (m):", res, "nodata:", nodata)

z = z.astype(np.float64)
mask = (z == nodata) | ~np.isfinite(z) | (z < -100)
z[mask] = np.nan

valid_mask = np.isfinite(z)
print("Valid (in-watershed) cells:", valid_mask.sum(), "/", valid_mask.size,
      "-> area (km^2):", valid_mask.sum() * res * res / 1e6)

# Fill small interior NaN gaps (holes fully surrounded by valid data) so a
# few stray no-data pixels inside the watershed don't break flow routing.
from scipy.ndimage import binary_fill_holes

filled_valid = binary_fill_holes(valid_mask)
interior_holes = filled_valid & ~valid_mask
if interior_holes.sum() > 0:
    print("Filling", interior_holes.sum(), "interior no-data cells with local mean")
    zz = z.copy()
    for _ in range(5):
        nanmask = ~np.isfinite(zz) & filled_valid
        if not nanmask.any():
            break
        # simple 3x3 nanmean fill, iterated
        padded = np.pad(zz, 1, mode="edge")
        for i, j in zip(*np.where(nanmask)):
            window = padded[i:i + 3, j:j + 3]
            with np.errstate(invalid="ignore"):
                m = np.nanmean(window)
            if np.isfinite(m):
                zz[i, j] = m
    z = zz
    valid_mask = filled_valid

np.savez(os.path.join(OUT_DIR, "dem_1m.npz"), z=z, valid=valid_mask, res=res)
print("Saved dem_1m.npz  (shape:", z.shape, ", res:", res, "m)")
print("NOTE: grid has", valid_mask.sum(), "valid nodes -- roughly",
      round(valid_mask.sum() / 11564, 1), "x more than the 15 m sweep's ~11,564 "
      "core nodes. Flow-routing cost per step scales at least linearly with "
      "node count (worse for the pure-Python depression finder), so budget "
      "HPC resources accordingly -- see README.md.")
