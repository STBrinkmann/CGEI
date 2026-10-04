# GAVI / lacunarity benchmark (works with CGEI 0.3.1 and >= 0.4.0).
#
# Usage: Rscript benchmarks/bench_gavi.R <data_dir> <out.rds> [sizes] [cores] [reps]
#   sizes  comma separated raster sizes, files gavi_<size>.tif (default 500,1000,2000)
#   cores  comma separated numbers of threads (default 1,4)
#   reps   repetitions, the median is reported (default 3)
#
# Times lacunarity(), the C++ focal step of gavi() (focal_sum) and the complete
# gavi() for a two-layer raster (binary greenspace + continuous layer).
suppressPackageStartupMessages({
  library(terra)
  library(CGEI)
})
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_file <- args[2]
sizes <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(500L, 1000L, 2000L)
cores_set <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else c(1L, 4L)
reps <- if (length(args) >= 5) as.integer(args[5]) else 3L

version <- as.character(utils::packageVersion("CGEI"))
new_api <- exists("raster_geometry", envir = asNamespace("CGEI"), inherits = FALSE)
geom <- function(r) if (new_api) CGEI:::raster_geometry(r) else raster::raster(terra::rast(r))
med_time <- function(expr_fun) stats::median(replicate(reps, system.time(expr_fun())[["elapsed"]]))

results <- list()
outputs <- list()
for (size in sizes) {
  x <- rast(file.path(data_dir, sprintf("gavi_%d.tif", size)))
  for (cores in cores_set) {
    lac <- NULL
    t_lac <- med_time(function() lac <<- lacunarity(x, cores = cores))
    xm <- values(x, mat = TRUE) * 1.0
    xg <- geom(x)
    lm <- as.matrix(lac[, c("i", "r", "Lac")])
    fs <- NULL
    t_focal <- med_time(function() fs <<- CGEI:::focal_sum(xg, xm, lm, TRUE, cores, FALSE))
    g <- NULL
    t_gavi <- med_time(function() {
      set.seed(1)
      g <<- gavi(x, lac, cores = cores)
    })
    results[[length(results) + 1]] <- data.frame(
      version = version, size = size, cores = cores, windows = paste(unique(lac$r), collapse = ","),
      lacunarity_s = t_lac, focal_s = t_focal, gavi_s = t_gavi)
    outputs[[paste(size, cores)]] <- list(lac = lac, focal = fs, gavi = values(g, mat = FALSE))
    cat(sprintf("%s %4d^2 cores=%d: lacunarity %8.3f s | focal step %8.3f s | gavi() %8.3f s\n",
                version, size, cores, t_lac, t_focal, t_gavi))
  }
}
saveRDS(list(timing = do.call(rbind, results), outputs = outputs), out_file)
