# VGVI benchmark (works with CGEI 0.3.1 and >= 0.4.0).
#
# Usage: Rscript benchmarks/bench_vgvi.R <data_dir> <out.rds> [n_obs] [reps] [radii] [cores]
#   data_dir  output of benchmarks/make_data.R
#   n_obs     number of observers per scene (random subset, default 1000)
#   reps      repetitions, the median is reported (default 3)
#   radii     comma separated max_distance values (default 100,200,300)
#   cores     comma separated numbers of threads (default 1,4)
#
# Times the C++ core (VGVI_cpp with prepared inputs) and the end-to-end R
# function vgvi() (mode = "exponential", m = 1, b = 3).
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(CGEI)
})
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_file <- args[2]
n_obs <- if (length(args) >= 3) as.integer(args[3]) else 1000L
reps <- if (length(args) >= 4) as.integer(args[4]) else 3L
radii <- if (length(args) >= 5) as.numeric(strsplit(args[5], ",")[[1]]) else c(100, 200, 300)
cores_set <- if (length(args) >= 6) as.integer(strsplit(args[6], ",")[[1]]) else c(1L, 4L)

version <- as.character(utils::packageVersion("CGEI"))
new_api <- exists("raster_geometry", envir = asNamespace("CGEI"), inherits = FALSE)
geom <- function(r) if (new_api) CGEI:::raster_geometry(r) else raster::raster(terra::rast(r))
med_time <- function(expr_fun) stats::median(replicate(reps, system.time(expr_fun())[["elapsed"]]))

results <- list()
values <- list()
for (scene in c("city", "open")) {
  dsm <- rast(file.path(data_dir, paste0(scene, "_dsm.tif")))
  dtm <- rast(file.path(data_dir, paste0(scene, "_dtm.tif")))
  gs <- rast(file.path(data_dir, paste0(scene, "_greenspace.tif")))
  obs <- st_read(file.path(data_dir, paste0("observers_", scene, ".gpkg")), quiet = TRUE)
  set.seed(1)
  obs <- obs[sort(sample.int(nrow(obs), min(n_obs, nrow(obs)))), ]
  xy <- st_coordinates(obs)
  h0 <- extract(dtm, xy)[, 1] + 1.7

  for (radius in radii) {
    # inputs prepared like in vgvi(): DSM / greenspace cropped to the observers' AOI
    aoi <- vect(st_buffer(st_as_sfc(st_bbox(obs)), radius + 2 * res(dsm)[1]))
    d <- crop(dsm, aoi, snap = "out")
    g <- crop(gs, aoi, snap = "out")
    dv <- values(d, mat = FALSE)
    gv <- as.numeric(values(g, mat = FALSE))
    dg <- geom(d)
    gg <- geom(g)
    c0 <- colFromX(d, xy[, 1])
    r0 <- rowFromY(d, xy[, 2])
    for (cores in cores_set) {
      v <- NULL
      t_core <- med_time(function() v <<- CGEI:::VGVI_cpp(dg, dv, gg, gv, c0, r0, h0, radius,
                                                          2L, 1, 3, cores, FALSE))
      t_e2e <- med_time(function() suppressMessages(
        vgvi(obs, dsm, dtm, gs, max_distance = radius, mode = "exponential", m = 1, b = 3, cores = cores)))
      results[[length(results) + 1]] <- data.frame(
        version = version, scene = scene, radius = radius, cores = cores, n_obs = nrow(obs),
        core_s = t_core, e2e_s = t_e2e, ms_per_obs = 1000 * t_core / nrow(obs),
        mean_vgvi = mean(v, na.rm = TRUE))
      values[[paste(scene, radius, cores)]] <- v
      cat(sprintf("%s %-4s r=%3d cores=%d: core %7.3f s (%.3f ms/obs), vgvi() %7.3f s\n",
                  version, scene, radius, cores, t_core, 1000 * t_core / nrow(obs), t_e2e))
    }
  }
}
saveRDS(list(timing = do.call(rbind, results), values = values), out_file)
