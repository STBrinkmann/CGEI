# Exactness at scale: compares the C++ results of CGEI (>= 0.4.0) on the
# benchmark data with the naive R reference implementation used by the tests.
#
# Usage: Rscript benchmarks/validate.R <data_dir> [n_obs] [max_distance]
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(CGEI)
})
source(file.path("tests", "testthat", "helper-reference.R"))
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
n_obs <- if (length(args) >= 2) as.integer(args[2]) else 15L
radius <- if (length(args) >= 3) as.numeric(args[3]) else 100

for (scene in c("city", "open")) {
  dsm <- rast(file.path(data_dir, paste0(scene, "_dsm.tif")))
  dtm <- rast(file.path(data_dir, paste0(scene, "_dtm.tif")))
  gs <- rast(file.path(data_dir, paste0(scene, "_greenspace.tif")))
  obs <- st_read(file.path(data_dir, paste0("observers_", scene, ".gpkg")), quiet = TRUE)
  set.seed(3)
  obs <- obs[sample.int(nrow(obs), n_obs), ]
  xy <- st_coordinates(obs)
  rows <- rowFromY(dsm, xy[, 2])
  cols <- colFromX(dsm, xy[, 1])
  h0 <- extract(dtm, xy)[, 1] + 1.7
  # window around the observers so the reference does not have to read the full raster
  r_cells <- round(radius / res(dsm)[1])
  r0 <- max(1, min(rows) - r_cells - 2); r1 <- min(nrow(dsm), max(rows) + r_cells + 2)
  c0 <- max(1, min(cols) - r_cells - 2); c1 <- min(ncol(dsm), max(cols) + r_cells + 2)
  win <- ext(xFromCol(dsm, c0) - 0.5, xFromCol(dsm, c1) + 0.5, yFromRow(dsm, r1) - 0.5, yFromRow(dsm, r0) + 0.5)
  dsm_w <- crop(dsm, win)
  gs_w <- crop(gs, win)
  dm <- as.matrix(dsm_w, wide = TRUE)
  gm <- as.matrix(gs_w, wide = TRUE)
  rr <- rowFromY(dsm_w, xy[, 2])
  cc <- colFromX(dsm_w, xy[, 1])

  t_ref <- system.time(ref <- sapply(seq_len(n_obs), function(k) {
    vs <- ref_viewshed(dm, rr[k], cc[k], h0[k], r_cells)
    ref_vgvi_from_rings(ref_rings(vs$visible, gm, rr[k], cc[k], res(dsm)[1]), radius, "exponential", 1, 3)
  }))[["elapsed"]]
  t_cpp <- system.time(v <- suppressMessages(
    vgvi(obs, dsm, dtm, gs, max_distance = radius, mode = "exponential", m = 1, b = 3)))[["elapsed"]]
  cat(sprintf("%s: %d observers, max_distance %g m: max |C++ - R reference| = %.2e (R reference %.1f s, vgvi() %.2f s)\n",
              scene, n_obs, radius, max(abs(v$VGVI - ref)), t_ref, t_cpp))
}
