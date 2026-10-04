# Benchmarks

Reproducible runtime comparisons of CGEI 0.3.1 (before the rewrite), 0.4.0
and 0.4.1. The folder is excluded from the package build (`.Rbuildignore`).

```sh
# 1. synthetic data (about 200 MB, written to benchmarks/data)
Rscript benchmarks/make_data.R benchmarks/data

# 2. install both versions into separate libraries, e.g.
git worktree add /tmp/cgei_old b7d8f34          # last commit of 0.3.1
R CMD INSTALL --library=lib_old /tmp/cgei_old
R CMD INSTALL --library=lib_new .

# 3. run every version in its own R process
#    (one after the other, nothing else running on the machine)
R_LIBS=lib_new Rscript benchmarks/bench_vgvi.R benchmarks/data benchmarks/results/vgvi_new.rds 1000 5
R_LIBS=lib_old Rscript benchmarks/bench_vgvi.R benchmarks/data benchmarks/results/vgvi_old.rds 1000 5
R_LIBS=lib_new Rscript benchmarks/bench_vvi.R benchmarks/data benchmarks/results/vvi_new.rds 500 3
R_LIBS=lib_old Rscript benchmarks/bench_vvi.R benchmarks/data benchmarks/results/vvi_old.rds 500 3
R_LIBS=lib_new Rscript benchmarks/bench_gavi.R benchmarks/data benchmarks/results/gavi_new.rds
# 0.3.1 needs minutes per call: one repetition, 1000 x 1000 with 4 threads only
R_LIBS=lib_old Rscript benchmarks/bench_gavi.R benchmarks/data benchmarks/results/gavi_old_500.rds 500 1,4 1
R_LIBS=lib_old Rscript benchmarks/bench_gavi.R benchmarks/data benchmarks/results/gavi_old_1000.rds 1000 4 1

# 4. tables
Rscript benchmarks/summarise.R benchmarks/results

# 5. exactness at scale: C++ results vs the naive R reference of the tests
R_LIBS=lib_new Rscript benchmarks/validate.R benchmarks/data 15 100
```

0.4.0 vs 0.4.1 (same scripts; `compare.R` prints the tables of RESULTS.md and
checks that the results are identical):

```sh
git worktree add /tmp/cgei_040 0.4.0
R CMD INSTALL --library=lib_040 /tmp/cgei_040
R CMD INSTALL --library=lib_041 .
for v in 040 041; do
  R_LIBS=lib_$v Rscript benchmarks/bench_vgvi.R benchmarks/data benchmarks/results/vgvi_$v.rds 1000 5 100,200,300 1,4
  R_LIBS=lib_$v Rscript benchmarks/bench_vgvi.R benchmarks/data benchmarks/results/vgvi_long_$v.rds 200 3 500,800 1,4
  R_LIBS=lib_$v Rscript benchmarks/bench_vvi.R benchmarks/data benchmarks/results/vvi_$v.rds 500 3
done
Rscript benchmarks/compare.R benchmarks/results 040 041 0.4.0 0.4.1
```

## Data (`make_data.R`, seeded)

* **city**: 2 x 2 km at 1 m resolution (4 million cells): undulating terrain,
  street grid with 100 m blocks, buildings of 6-30 m (a few towers up to 60 m),
  parks and street trees; greenspace = tree crowns and parks (19 % green).
  Observers every 10 m on the street centre lines (3,920 points).
* **open**: the same terrain without buildings and with scattered trees
  (worst case for the early termination of lines of sight); observers on a
  20 m grid (4,900 points).
* **gavi_500 / 1000 / 2000**: two-layer rasters (clustered binary greenspace
  and a continuous layer) for `lacunarity()` and `gavi()` with the default box
  sizes (up to half of the raster size).

The results of the last run are in [RESULTS.md](RESULTS.md).
