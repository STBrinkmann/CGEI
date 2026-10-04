# Stand-alone C++ tests (sanitizers)

The performance-critical code of CGEI lives in plain C++ headers without any
R / Rcpp dependency:

| header | used by |
|---|---|
| `src/los_geometry.h` | Bresenham lines of sight (`LoS_reference()`) |
| `src/viewshed_engine.h` | viewshed sweep of `vgvi()`, `vvi()`, `viewshed_list()`; potential viewshed and per-cell counts of `vvi()` |
| `src/boxfilter.h` | box sums / counts / max / min of `gavi()` and `lacunarity()` |
| `src/natural_breaks.h` | Jenks / Fisher natural breaks of `gavi()` |

`engine_tests.cpp` compiles exactly these headers, runs them multi-threaded
(OpenMP) on random data and compares every result with a naive single-threaded
re-implementation. The batched viewshed sweep is checked with and without early
termination, with double and float DSM storage and with batch sizes 1, 3 and
16 (masks must come out in raster order and cleared). This includes the per-cell
visibility counts of `vvi(mode = "cumulative" / "viewshed")`, which several
threads accumulate concurrently (atomic increments). It also contains a structural port of the
original (CGEI 0.3.1) VGVI/VVI viewshed loop to check its OpenMP structure.

```sh
make release   # -O2, plain run
make asan      # AddressSanitizer + UndefinedBehaviorSanitizer (gcc + libgomp)
make tsan      # ThreadSanitizer (clang + LLVM OpenMP runtime + Archer)
```

Requirements (Ubuntu): `g++`, `clang`, `libomp-dev`, `libclang-rt-dev`.

ThreadSanitizer needs the Archer tool of the LLVM OpenMP runtime
(`OMP_TOOL_LIBRARIES=.../libarcher.so`) and
`TSAN_OPTIONS=ignore_noninstrumented_modules=1`, otherwise it reports false
positives inside the (uninstrumented) OpenMP runtime. The Makefile sets both.
With this setup a deliberate race is reported and race-free code is clean.

## Running the R test suite under ASan / UBSan

```sh
cat > Makevars.asan <<'EOF'
CXXFLAGS = -g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=undefined
CXX17FLAGS = $(CXXFLAGS)
LDFLAGS = -fsanitize=address,undefined
EOF
R_MAKEVARS_USER=$PWD/Makevars.asan R CMD INSTALL --preclean --no-test-load --library=lib_asan .
cd tests
LD_PRELOAD="$(gcc -print-file-name=libasan.so) $(gcc -print-file-name=libubsan.so)" \
ASAN_OPTIONS=detect_leaks=0:abort_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
R_LIBS=../lib_asan Rscript -e 'testthat::test_dir("testthat")'
```

## Last results (2026-10-04, CGEI 0.4.1, gcc 13.3 / clang 18.1, 4 threads)

| run | result |
|---|---|
| `make release` | 767 checks, 0 failures |
| `make asan` | 767 checks, 0 failures, no sanitizer reports |
| `make tsan` | 767 checks, 0 failures, no data races (new engines and old loop) |
| R test suite with ASan/UBSan-instrumented package | 504 expectations pass, no reports |

The fast Fisher natural breaks (`FisherDC`, O(k n log n)) give the same
partition as the exact port of classInt's Fortran routine in all cases with
continuous data; with many repeated values (near-ties) they may pick a
different partition, which is never worse (its within-class sum of squares is
equal or lower, the difference being rounding noise of the Fortran routine).
