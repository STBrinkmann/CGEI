# Benchmark results

## CGEI 0.4.0 vs 0.4.1

Run on 2026-10-04 with the scripts in this folder (see [README.md](README.md)).

* Machine: cloud VM (Firecracker), Intel Xeon (Cascade Lake) @ 2.8 GHz, 4 vCPUs,
  15 GB RAM, Ubuntu 24.04. A different (and somewhat slower) VM than the one
  of the 0.3.1 vs 0.4.0 comparison below, so absolute times differ.
* R 4.3.3, gcc 13.3 (`-O2`, OpenMP), terra 1.7-65, Rcpp 1.0.12
* 0.4.0 = release 0.4.0 (commit `33a6632`), 0.4.1 = this branch
* All benchmarks ran one after the other with nothing else running on the
  machine; both versions in separate R processes. Medians of 5 (VGVI 100-300
  m) or 3 (VGVI 500-800 m, VVI) runs.
* Same data as below (`make_data.R`): 1000 (100-300 m) / 200 (500-800 m)
  randomly chosen observers of the city and open-terrain scenes,
  `mode = "exponential"` (m = 1, b = 3), 1 m DSM (Float32 GeoTIFF).
* All VGVI and VVI results are identical (`identical()`) between both
  versions.
* The complete calls contain about 0.2-0.4 s of cropping and reading the
  4-million-cell rasters with `terra`, which varies by about +-20 % between
  repeated runs (e.g. `vgvi()` city 100 m with 4 threads: 0.37-0.46 s for
  both versions in repeated runs, while the C++ core takes 0.08 s). For small
  radii the complete calls are dominated by this part.

### Overview

| function | 1 thread | 4 threads |
|---|---|---|
| `vgvi()`, C++ core, `max_distance` 100-300 m | 1.2-2.2x | 1.0-1.6x |
| `vgvi()`, C++ core, `max_distance` 500-800 m | 1.7-2.2x | 1.4-1.8x |
| `vgvi()`, complete call | 1.0-2.2x | 1.0-1.6x |
| `vvi()`, complete call (all modes) | 0.9-1.3x | 0.9-1.1x |

The gain grows with `max_distance`: the larger the neighbourhood of an
observer, the more the old sweep was limited by memory accesses. With 4
threads the gain is smaller, because the old code profited more from the
additional per-core caches (it scaled super-linearly, up to 4.7x on 4
threads). `vvi()` changes little: its complete calls are dominated by reading
the rasters.

What the time was spent on (C++ core; `perf` cpu-clock sampling and timing
experiments on the same data):

* The sweep is limited by memory accesses: on an artificial surface on which
  every visibility test is perfectly predictable, a step of a line of sight
  costs 5 cycles for r = 100 cells and 12-15 cycles for r = 300-600 cells; it
  drops back to about 6 cycles if either the line-of-sight table or the DSM is
  forced into the L1 cache. The division and branch mispredictions are minor.
* In `vgvi()`, 30-50 % of the time went into the bookkeeping of the visible
  cells (byte mask, list of touched cells, ring and greenspace look-ups with
  random access); 0.4.1 replaces it by one bit operation per visible cell and
  a sequential pass over the mask per observer.
* The lines of sight process 2.0-2.1 steps per cell of the circle (the
  shared-prefix reuse saves only 6-18 % of the steps, as neighbouring lines
  split up and then run over the same cells again). This is inherent to the
  "visible from any line of sight" definition and was not changed.
* Tried and discarded (exact, but not faster on this machine): skipping
  blocks of 8 x 8 cells that cannot be visible (reduced the evaluated steps by
  2-3x in the city, but the additional bound checks cost as much), a
  branch-free inner loop, SIMD over 4 observers in lock step, bit-encoded
  lines of sight with distances computed on the fly, a transposed DSM copy for
  north-south lines, and sweeping the 8 symmetric octants of a base line
  together.

### VGVI: C++ core (VGVI_cpp), milliseconds per observer (1000 observers)

| scene | max_distance | threads | 0.4.0 | 0.4.1 | speed-up |
|---|---|---|---|---|---|
| city | 100 | 1 | 0.349 | 0.277 | **1.3x** |
| city | 100 | 4 | 0.084 | 0.084 | **1.0x** |
| city | 200 | 1 | 0.912 | 0.672 | **1.4x** |
| city | 200 | 4 | 0.242 | 0.211 | **1.1x** |
| city | 300 | 1 | 2.162 | 1.270 | **1.7x** |
| city | 300 | 4 | 0.465 | 0.351 | **1.3x** |
| open | 100 | 1 | 0.539 | 0.463 | **1.2x** |
| open | 100 | 4 | 0.141 | 0.114 | **1.2x** |
| open | 200 | 1 | 2.177 | 1.114 | **2.0x** |
| open | 200 | 4 | 0.464 | 0.353 | **1.3x** |
| open | 300 | 1 | 4.672 | 2.096 | **2.2x** |
| open | 300 | 4 | 0.993 | 0.621 | **1.6x** |

### VGVI: complete vgvi() call for 1000 observers, seconds

| scene | max_distance | threads | 0.4.0 | 0.4.1 | speed-up |
|---|---|---|---|---|---|
| city | 100 | 1 | 0.60 | 0.57 | **1.0x** |
| city | 100 | 4 | 0.37 | 0.46 | **0.8x** |
| city | 200 | 1 | 1.17 | 1.00 | **1.2x** |
| city | 200 | 4 | 0.54 | 0.55 | **1.0x** |
| city | 300 | 1 | 2.41 | 1.53 | **1.6x** |
| city | 300 | 4 | 0.96 | 0.64 | **1.5x** |
| open | 100 | 1 | 0.81 | 0.63 | **1.3x** |
| open | 100 | 4 | 0.42 | 0.38 | **1.1x** |
| open | 200 | 1 | 2.63 | 1.49 | **1.8x** |
| open | 200 | 4 | 0.74 | 0.70 | **1.1x** |
| open | 300 | 1 | 5.06 | 2.31 | **2.2x** |
| open | 300 | 4 | 1.35 | 0.88 | **1.5x** |

VGVI values identical in 12 of 12 configurations.

### VGVI: C++ core (VGVI_cpp), milliseconds per observer (200 observers)

| scene | max_distance | threads | 0.4.0 | 0.4.1 | speed-up |
|---|---|---|---|---|---|
| city | 500 | 1 | 6.585 | 3.045 | **2.2x** |
| city | 500 | 4 | 2.160 | 1.200 | **1.8x** |
| city | 800 | 1 | 11.330 | 6.640 | **1.7x** |
| city | 800 | 4 | 4.395 | 3.135 | **1.4x** |
| open | 500 | 1 | 10.540 | 5.615 | **1.9x** |
| open | 500 | 4 | 3.115 | 1.835 | **1.7x** |
| open | 800 | 1 | 21.645 | 11.110 | **1.9x** |
| open | 800 | 4 | 6.535 | 4.155 | **1.6x** |

### VGVI: complete vgvi() call for 200 observers, seconds

| scene | max_distance | threads | 0.4.0 | 0.4.1 | speed-up |
|---|---|---|---|---|---|
| city | 500 | 1 | 1.69 | 0.96 | **1.8x** |
| city | 500 | 4 | 0.76 | 0.52 | **1.4x** |
| city | 800 | 1 | 2.62 | 1.60 | **1.6x** |
| city | 800 | 4 | 1.31 | 0.91 | **1.4x** |
| open | 500 | 1 | 2.61 | 1.45 | **1.8x** |
| open | 500 | 4 | 1.01 | 0.69 | **1.5x** |
| open | 800 | 1 | 4.69 | 2.72 | **1.7x** |
| open | 800 | 4 | 1.64 | 1.05 | **1.6x** |

VGVI values identical in 8 of 8 configurations.

### VVI (city, 500 observers), seconds

| max_distance | threads | vvi() 0.4.0 | vvi() 0.4.1 | speed-up | cumulative 0.4.0 | cumulative 0.4.1 | speed-up  | viewshed 0.4.0 | viewshed 0.4.1 |  speed-up |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 1 | 0.40 | 0.41 | **1.0x** | 0.42 | 0.45 | **0.9x** | 0.69 | 0.75 | **0.9x** |
| 100 | 4 | 0.21 | 0.22 | **1.0x** | 0.24 | 0.22 | **1.1x** | 0.48 | 0.47 | **1.0x** |
| 200 | 1 | 0.66 | 0.53 | **1.3x** | 0.78 | 0.64 | **1.2x** | 1.32 | 1.06 | **1.2x** |
| 200 | 4 | 0.31 | 0.33 | **0.9x** | 0.38 | 0.34 | **1.1x** | 0.90 | 0.86 | **1.0x** |

VVI results (mean VVI, cumulative VVI, sum of n_views) identical in 4 of 4 configurations.

## CGEI 0.3.1 vs 0.4.0

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

### Overview

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

#### VGVI: C++ core (VGVI_cpp), milliseconds per observer

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

#### VGVI: complete vgvi() call for 1000 observers, seconds

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

#### VGVI: mean value of the index (results change because of the bug fixes)

| scene | max_distance | 0.3.1 | 0.4.0 |
|---|---|---|---|
| city | 100 | 0.1195 | 0.1364 |
| city | 200 | 0.1359 | 0.1575 |
| city | 300 | 0.1436 | 0.1744 |
| open | 100 | 0.3511 | 0.3572 |
| open | 200 | 0.3570 | 0.3639 |
| open | 300 | 0.3650 | 0.3731 |

#### GAVI: seconds (two-layer raster, default lacunarity box sizes)

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

#### VVI (city, 500 observers), seconds

| max_distance | threads | vvi() 0.3.1 | vvi() 0.4.0 | speed-up | cumulative 0.3.1 | cumulative 0.4.0 |  speed-up | viewshed 0.3.1 | viewshed 0.4.0 | speed-up   |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 1 | 7.29 | 0.22 | **33x** | 8.38 | 0.25 | **33x** | 8.26 | 0.43 | **19x** |
| 100 | 4 | 3.26 | 0.15 | **22x** | 4.39 | 0.16 | **27x** | 3.81 | 0.27 | **14x** |
| 200 | 1 | 37.28 | 0.36 | **102x** | 37.61 | 0.42 | **89x** | 31.47 | 0.73 | **43x** |
| 200 | 4 | 13.62 | 0.19 | **73x** | 15.01 | 0.22 | **69x** | 13.50 | 0.52 | **26x** |

#### VVI: mean values (results change slightly because of the bug fixes)

| max_distance | mean VVI 0.3.1 | mean VVI 0.4.0 | cumulative VVI 0.3.1 | cumulative VVI 0.4.0 |
|---|---|---|---|---|
| 100 | 0.2326 | 0.2323 | 0.5940 | 0.5938 |
| 200 | 0.0785 | 0.0784 | 0.5084 | 0.5083 |

