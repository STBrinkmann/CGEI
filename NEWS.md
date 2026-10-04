# CGEI 0.4.0

This release fixes several bugs that affected the results of `vgvi()`,
`vvi()`, `viewshed_list()`, `lacunarity()` and `gavi()`, and makes them much
faster (see `benchmarks/RESULTS.md`). **VGVI values change** (substantially for
rasters with a resolution above 1 m and for observers near the right border
of the cropped DSM), VVI values change slightly.

## Bug fixes

### `vgvi()`

-   Greenspace values and distances were taken from the wrong cell: the viewshed
    stored 1-based cell numbers that were then used as 0-based indices. Every
    visible cell was evaluated one column to the east, the observer cell got a
    distance of one cell, and visible cells in the last raster column wrapped to
    the first column of the next row (distance = raster width), which pushed
    VGVI towards 0. The test value of the package changed from 0.004 / 0.001 to
    0.271 / 0.5 (random greenspace with 50 % green cells).
-   Distance rings (1 map unit) without any visible cell counted as 0 % green
    but kept their decay weight. They are now ignored, i.e. VGVI = 1 if all
    visible cells are green, independent of the raster resolution (a 100 %
    green view gave 0.91 at 5 m resolution before).
-   Fractional greenspace values were truncated to integers, and a viewshed
    consisting of a single distance ring used integer division.

### Line of sight (`vgvi()`, `vvi()`, `viewshed_list()`)

-   An obstacle in the first cell of a line of sight was ignored whenever that
    line shared exactly this first cell with the previous line (`k_i > 1`
    instead of `k_i > 0`). An observer enclosed by a 10 m wall saw 249 instead
    of 9 cells.
-   Cells with NA heights skipped the horizon book-keeping, so later lines could
    reuse a stale horizon of an unrelated line.
-   Rasters with at most 2 * `max_distance` / resolution columns (e.g. a single
    observer) could let cells from the opposite raster border enter the
    viewshed (column wrap-around).
-   The DSM is cropped with `snap = "out"` and an extra margin, so observers at
    the border of the area of interest keep their full circle.
-   `vvi(by_row = TRUE)` failed if observer points were removed (outside of the
    DSM / DTM).

### `lacunarity()` and `gavi()`

-   `lacunarity()`: a box whose newly added rim contained only NA became NA and
    lost the mass of its valid core; results depended on the order of `r_vec`;
    box sizes larger than the raster read outside of the raster (they are now
    dropped with a warning).
-   `gavi()`: layers with fewer than 9 distinct values produced a broken
    reclassification matrix; sampling warned for rasters with NA cells.

## Performance

Speed-up compared with CGEI 0.3.1 on generated test data (1 m city and
open-terrain scenes of 2 x 2 km, 500 x 500 to 1000 x 1000 GAVI rasters; 1 and
4 threads; details and reproduction in `benchmarks/RESULTS.md`):

| function | speed-up |
|---|---|
| `vgvi()`, C++ core | 6-10x |
| `vgvi()`, complete call (1000 observers) | 6-9x |
| `vvi()` | 22-102x |
| `vvi(mode = "cumulative")` / `vvi(mode = "viewshed")` | 27-89x / 14-43x |
| `lacunarity()` | 156-516x |
| `gavi()` | 345-1557x |

-   New C++ viewshed engine: precomputed line-of-sight geometry and decay
    weights, no memory allocation per observer, exact early termination of
    lines of sight that cannot reveal anything anymore (block maxima of the
    DSM), observers processed in a cache-friendly (Morton) order. No forced
    garbage collection in `vgvi()`, `vvi()` and `viewshed_list()`.
-   `vvi()`: the default mode only counts cells; the cumulative and viewshed
    modes accumulate per-cell counts in C++ instead of returning the cells of
    every observer to R (memory: two integer rasters instead of one vector per
    observer).
-   `gavi()` focal means and `lacunarity()` box masses use separable sliding
    windows: O(1) per cell instead of O(window size^2). Results are
    bit-identical to 0.3.1 for integer and single-precision (GeoTIFF float)
    rasters.
-   The Jenks / Fisher natural breaks of `gavi()` are computed in C++ (identical
    breaks to `classInt`, Fisher's optimum in O(k n log n) instead of
    O(k n^2)), and the reclassification runs on the values in memory instead
    of `terra` round trips.

## Multi-threading

-   No R API calls or Rcpp objects inside OpenMP regions anymore (an
    out-of-range access in a worker thread could crash R), thread-safe progress
    bar, `num_threads()` instead of changing the global OpenMP setting, and
    `cores` is validated (a warning is given if CGEI was built without OpenMP).
    Results are identical for any number of threads.

## Other changes

-   The `raster` and `classInt` packages are no longer required.
-   `gavi()` draws the cells for the natural breaks (at most 50,000) with
    `sample.int()` instead of `terra::spatSample()`. Rasters with at most
    50,000 valid cells use all cells, as before; for larger rasters,
    `set.seed()` makes the result reproducible.
-   New test suite with shipped test data, an independent R reference
    implementation and tests for R/C++ index offsets and multi-threading;
    stand-alone sanitizer tests in `dev/cpp-tests`, benchmarks in `benchmarks/`.

# CGEI 0.3.1

-   Fixed the "subscript out of bounds" bug in the `vgvi.cpp` and `vvi.cpp` functions.

# CGEI 0.3.0

## New Features

-   `gavi()` function has been added. This function calculates the Greenspace Availability IndexIndex (GAVI) after a `lacunarity()` analysis.

## Docs

-   The C++ `lacunarity()` function now uses the R4 raster class.
-   The parameters `cores` and `progress` were adjusted in a uniform style in the existing functions.
-   Tests for the `gavi()` function have been added.


# CGEI 0.2.1

## Bug fixes

-   `vvi(mode = c("VVI", "cumulative"), by_row = TRUE)` now works as expected. ([#17](https://github.com/STBrinkmann/CGEI/issues/17))
-   `vvi()` manual has been adjusted. ([#23](https://github.com/STBrinkmann/CGEI/issues/23))

# CGEI 0.2.0

## New Features

-   `sf_interpolate_IDW()` function has been added. This is a highly efficient IDW method for converting an `sf` object to a `SpatRast`.

# CGEI 0.1.0

This is the first release of CGEI.

## New Features

-   `vvi()` function now accepts `mode = "cumulative"` and `by_row` arguments. ([#17](https://github.com/STBrinkmann/CGEI/issues/17), @hansvancalster)
