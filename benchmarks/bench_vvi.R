# VVI benchmark (works with CGEI 0.3.1 and >= 0.4.0).
# Usage: Rscript benchmarks/bench_vvi.R <data_dir> <out.rds> [n_obs] [reps]
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
version <- as.character(utils::packageVersion("CGEI"))
med_time <- function(expr_fun) stats::median(replicate(reps, system.time(expr_fun())[["elapsed"]]))

dsm <- rast(file.path(data_dir, "city_dsm.tif"))
dtm <- rast(file.path(data_dir, "city_dtm.tif"))
obs <- st_read(file.path(data_dir, "observers_city.gpkg"), quiet = TRUE)
set.seed(1)
obs <- obs[sort(sample.int(nrow(obs), n_obs)), ]
results <- list()
for (radius in c(100, 200)) {
  for (cores in c(1L, 4L)) {
    v <- NULL
    t_vvi <- med_time(function() v <<- suppressMessages(vvi(obs, dsm, dtm, max_distance = radius, cores = cores)))
    cvvi <- NULL
    t_cum <- med_time(function() cvvi <<- suppressMessages(vvi(obs, dsm, dtm, max_distance = radius,
                                                               mode = "cumulative", cores = cores)))
    vs <- NULL
    t_vs <- med_time(function() vs <<- suppressMessages(vvi(obs, dsm, dtm, max_distance = radius,
                                                            mode = "viewshed", cores = cores)))
    results[[length(results) + 1]] <- data.frame(version = version, radius = radius, cores = cores,
                                                 n_obs = n_obs, vvi_s = t_vvi, cumulative_s = t_cum,
                                                 viewshed_s = t_vs, mean_vvi = mean(v$VVI), cvvi = cvvi,
                                                 sum_n_views = sum(terra::values(vs$n_views), na.rm = TRUE))
    cat(sprintf("%s r=%d cores=%d: vvi() %.3f s, mode = \"cumulative\" %.3f s, mode = \"viewshed\" %.3f s\n",
                version, radius, cores, t_vvi, t_cum, t_vs))
  }
}
saveRDS(list(timing = do.call(rbind, results)), out_file)
