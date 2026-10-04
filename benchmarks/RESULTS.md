# Benchmark results: CGEI 0.3.1 vs 0.4.0

Run on 2026-10-04 with the scripts in this folder (see [README.md](README.md)).

* Machine: cloud VM (Firecracker), Intel Xeon @ 2.1 GHz, 4 vCPUs, 15 GB RAM, Ubuntu 24.04
* R 4.3.3, gcc 13.3 (`-O2`, OpenMP), terra 1.7-65, Rcpp 1.0.12
* 0.3.1 = commit `b7d8f34` (main), 0.4.0 = this branch
* All benchmarks ran one after the other with nothing else running on the
  machine; old and new version in separate R processes.
* Every number is the median of 5 (VGVI) or 3 (VVI, GAVI) runs; the slow 0.3.1
  GAVI cases ran once.
* Repeated runs agree to within a few percent single-threaded and to within
  about 10 % with 4 threads.
* VGVI: 1000 randomly chosen observers of each scene, `mode = "exponential"`
  (m = 1, b = 3), 1 m DSM. "C++ core" = `VGVI_cpp()` on prepared inputs,
  "complete call" = `vgvi()` including cropping, extraction etc.
* VVI: 500 randomly chosen observers of the city scene.
* GAVI: `lacunarity()` with the default box sizes, the focal step of `gavi()`
  (`focal_sum()`) and the complete `gavi()` (focal step and Jenks
  reclassification) for a two-layer raster (clustered binary greenspace and a
  continuous layer). 0.3.1 was not run for 2000 x 2000 (its focal step alone
  would take hours).
* Exactness at scale (`validate.R`): for 15 observers per scene, the VGVI of
  0.4.0 equals the naive R reference implementation of the tests to within
  4.4e-16 (max_distance 100 m and 200 m).

## Overview

Speed-up of 0.4.0 over 0.3.1 (ranges over the scenes, distances and raster
sizes of the tables below):

| function | 1 thread | 4 threads |
|---|---|---|
| `vgvi()`, C++ core | 7-10x | 6-9x |
| `vgvi()`, complete call (1000 observers) | 7-9x | 6-7x |
| `vvi()` | 33-102x | 22-73x |
| `vvi(mode = "cumulative")` | 33-89x | 27-69x |
| `vvi(mode = "viewshed")` | 19-43x | 14-26x |
| `lacunarity()` | 250x | 156-516x |
| `gavi()`, focal step | 4885x | 2383-10147x |
| `gavi()`, complete call | 703x | 345-1557x |

* VGVI: the work per observer is the visibility sweep over the cells within
  `max_distance` (0.4.0, single-threaded: 6-10 ns per cell of the circle,
  e.g. 0.19-0.30 ms per observer at 100 m).
  0.4.0 removes the per-observer allocations and the per-cell bookkeeping of
  0.3.1, precomputes the line-of-sight geometry and decay weights, ends lines
  of sight early when no further cell can be visible, and processes the
  observers in a cache-friendly order. 4 threads are 3.4-4.0 times faster
  than 1 thread. The complete call additionally contains about 0.15-0.2 s of
  cropping and reading the 4-million-cell rasters with `terra`.
* VVI: the default mode only counts cells; the cumulative and viewshed modes
  accumulate per-cell counts in C++ instead of returning the cells of every
  observer to R (for 500 observers at 200 m that were 62.8 million cell
  numbers, i.e. several hundred MB of R vectors).
* GAVI: the speed-up grows with the raster size, because 0.3.1 computed every
  focal mean and every lacunarity box in O(window size^2) with windows up to
  half the raster size. 0.3.1 would need hours for the 2000 x 2000 raster.
* Results: lacunarity values are identical to 0.3.1; focal means are identical
  for the binary layer and equal to within 5.6e-16 (relative) for the
  continuous layer. VGVI and VVI change because of the bug fixes (see the
  mean values below and `NEWS.md`).

### VGVI: C++ core (VGVI_cpp), milliseconds per observer

| scene | max_distance | threads | 0.3.1 | 0.4.0 | speed-up |
|---|---|---|---|---|---|
| city | 100 | 1 | 1.697 | 0.187 | **9x** |
| city | 100 | 4 | 0.408 | 0.055 | **7x** |
| city | 200 | 1 | 4.402 | 0.506 | **9x** |
| city | 200 | 4 | 1.083 | 0.137 | **8x** |
| city | 300 | 1 | 8.447 | 0.868 | **10x** |
| city | 300 | 4 | 1.963 | 0.226 | **9x** |
| open | 100 | 1 | 2.189 | 0.295 | **7x** |
| open | 100 | 4 | 0.557 | 0.085 | **7x** |
| open | 200 | 1 | 5.899 | 0.903 | **7x** |
| open | 200 | 4 | 1.420 | 0.226 | **6x** |
| open | 300 | 1 | 10.828 | 1.522 | **7x** |
| open | 300 | 4 | 2.588 | 0.396 | **7x** |

### VGVI: complete vgvi() call for 1000 observers, seconds

| scene | max_distance | threads | 0.3.1 | 0.4.0 | speed-up |
|---|---|---|---|---|---|
| city | 100 | 1 | 2.66 | 0.38 | **7x** |
| city | 100 | 4 | 1.31 | 0.23 | **6x** |
| city | 200 | 1 | 5.17 | 0.72 | **7x** |
| city | 200 | 4 | 1.93 | 0.32 | **6x** |
| city | 300 | 1 | 9.39 | 1.06 | **9x** |
| city | 300 | 4 | 2.87 | 0.42 | **7x** |
| open | 100 | 1 | 3.08 | 0.46 | **7x** |
| open | 100 | 4 | 1.38 | 0.25 | **6x** |
| open | 200 | 1 | 6.87 | 1.04 | **7x** |
| open | 200 | 4 | 2.36 | 0.40 | **6x** |
| open | 300 | 1 | 11.86 | 1.65 | **7x** |
| open | 300 | 4 | 3.41 | 0.62 | **6x** |

### VGVI: mean value of the index (results change because of the bug fixes)

| scene | max_distance | 0.3.1 | 0.4.0 |
|---|---|---|---|
| city | 100 | 0.1195 | 0.1364 |
| city | 200 | 0.1359 | 0.1575 |
| city | 300 | 0.1436 | 0.1744 |
| open | 100 | 0.3511 | 0.3572 |
| open | 200 | 0.3570 | 0.3639 |
| open | 300 | 0.3650 | 0.3731 |

### GAVI: seconds (two-layer raster, default lacunarity box sizes)

| raster | threads | box sizes | lacunarity 0.3.1 | lacunarity 0.4.0 | speed-up | focal step 0.3.1 | focal step 0.4.0 | speed-up  | gavi() 0.3.1 | gavi() 0.4.0 |  speed-up |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 500 x 500 | 1 | 3,5,9,17,33,65,129,251 | 32.53 | 0.13 | **250x** | 102.59 | 0.021 | **4885x** | 136.34 | 0.19 | **703x** |
| 500 x 500 | 4 | 3,5,9,17,33,65,129,251 | 9.06 | 0.06 | **156x** | 28.59 | 0.012 | **2383x** | 60.68 | 0.18 | **345x** |
| 1000 x 1000 | 1 | 3,5,9,17,33,65,129,257,501 | - | 0.49 | - | - | 0.103 | - | - | 0.36 | - |
| 1000 x 1000 | 4 | 3,5,9,17,33,65,129,257,501 | 129.00 | 0.25 | **516x** | 456.61 | 0.045 | **10147x** | 487.31 | 0.31 | **1557x** |
| 2000 x 2000 | 1 | 3,5,9,17,33,65,129,257,513,1001 | - | 2.62 | - | - | 0.536 | - | - | 1.26 | - |
| 2000 x 2000 | 4 | 3,5,9,17,33,65,129,257,513,1001 | - | 1.15 | - | - | 0.198 | - | - | 0.90 | - |

- 500 (threads 1): lacunarity identical: TRUE; focal step: layer 1 identical: TRUE, max. relative difference 0.0e+00
- 500 (threads 4): lacunarity identical: TRUE; focal step: layer 1 identical: TRUE, max. relative difference 0.0e+00
- 1000 (threads 4): lacunarity identical: TRUE; focal step: layer 1 identical: TRUE, max. relative difference 5.6e-16

### VVI (city, 500 observers), seconds

| max_distance | threads | vvi() 0.3.1 | vvi() 0.4.0 | speed-up | cumulative 0.3.1 | cumulative 0.4.0 |  speed-up | viewshed 0.3.1 | viewshed 0.4.0 | speed-up   |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 1 | 7.29 | 0.22 | **33x** | 8.38 | 0.25 | **33x** | 8.26 | 0.43 | **19x** |
| 100 | 4 | 3.26 | 0.15 | **22x** | 4.39 | 0.16 | **27x** | 3.81 | 0.27 | **14x** |
| 200 | 1 | 37.28 | 0.36 | **102x** | 37.61 | 0.42 | **89x** | 31.47 | 0.73 | **43x** |
| 200 | 4 | 13.62 | 0.19 | **73x** | 15.01 | 0.22 | **69x** | 13.50 | 0.52 | **26x** |

### VVI: mean values (results change slightly because of the bug fixes)

| max_distance | mean VVI 0.3.1 | mean VVI 0.4.0 | cumulative VVI 0.3.1 | cumulative VVI 0.4.0 |
|---|---|---|---|---|
| 100 | 0.2326 | 0.2323 | 0.5940 | 0.5938 |
| 200 | 0.0785 | 0.0784 | 0.5084 | 0.5083 |

