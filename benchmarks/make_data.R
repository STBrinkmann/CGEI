# Synthetic, reproducible benchmark data for CGEI.
#
# Usage: Rscript benchmarks/make_data.R [output_dir]   (default: benchmarks/data)
#
# Creates
#   city_{dsm,dtm,greenspace}.tif  1 m "city" of 2 x 2 km: gently undulating
#                                  terrain, street grid, building blocks
#                                  (6-30 m, a few towers up to 60 m), parks and
#                                  street trees (tree crowns are greenspace).
#   open_{dsm,dtm,greenspace}.tif  same terrain without buildings, sparse trees
#                                  (worst case for early termination).
#   observers_city.gpkg            points every 10 m along the street centre lines.
#   observers_open.gpkg            regular 20 m grid of points.
#   gavi_{500,1000,2000}.tif       two-layer rasters for lacunarity() / gavi():
#                                  clustered binary greenspace + continuous layer.
suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else file.path("benchmarks", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
crs <- "EPSG:25832"

set.seed(2024)
n <- 2000                      # cells per side (1 m resolution)
xmin <- 500000; ymax <- 5400000

# Smooth random field: sum of a few random cosine waves
smooth_field <- function(n, k = 6, wl = c(150, 800), amp = 1) {
  x <- seq_len(n)
  f <- matrix(0, n, n)
  for (i in seq_len(k)) {
    th <- runif(1, 0, 2 * pi); l <- runif(1, wl[1], wl[2]); ph <- runif(1, 0, 2 * pi)
    f <- f + outer(x, x, function(r, c) cos(2 * pi * (r * sin(th) + c * cos(th)) / l + ph))
  }
  amp * f / k
}

dtm <- 300 + smooth_field(n, amp = 12)

# Street grid: blocks of 100 m, streets 16 m wide
block <- 100; street <- 16
is_street_row <- (seq_len(n) - 1) %% block < street
is_street <- outer(is_street_row, is_street_row, `|`)

buildings <- matrix(0, n, n)
park <- matrix(FALSE, n, n)
nb <- n / block
for (bi in seq_len(nb)) {
  for (bj in seq_len(nb)) {
    r0 <- (bi - 1) * block + street + 1; r1 <- bi * block
    c0 <- (bj - 1) * block + street + 1; c1 <- bj * block
    if (runif(1) < 0.18) {            # park
      park[r0:r1, c0:c1] <- TRUE
      next
    }
    # 3-6 buildings per block
    for (k in seq_len(sample(3:6, 1))) {
      h <- if (runif(1) < 0.04) runif(1, 35, 60) else runif(1, 6, 30)
      w <- sample(15:40, 2)
      rr <- sample(r0:(r1 - w[1]), 1); cc <- sample(c0:(c1 - w[2]), 1)
      buildings[rr:(rr + w[1] - 1), cc:(cc + w[2] - 1)] <- h
    }
  }
}

# Tree crowns: many in parks, some along streets
add_trees <- function(n_trees, where, canopy) {
  idx <- which(where)
  centres <- idx[sample.int(length(idx), n_trees)]
  rows <- (centres - 1) %% n + 1; cols <- (centres - 1) %/% n + 1
  for (k in seq_along(centres)) {
    rad <- runif(1, 2, 6); h <- runif(1, 8, 20)
    rr <- max(1, floor(rows[k] - rad)):min(n, ceiling(rows[k] + rad))
    cc <- max(1, floor(cols[k] - rad)):min(n, ceiling(cols[k] + rad))
    d2 <- outer((rr - rows[k])^2, (cc - cols[k])^2, `+`)
    crown <- h * sqrt(pmax(0, 1 - d2 / rad^2))
    canopy[rr, cc] <- pmax(canopy[rr, cc], crown)
  }
  canopy
}
canopy <- matrix(0, n, n)
canopy <- add_trees(6000, park, canopy)
canopy <- add_trees(3000, is_street & !outer(abs((seq_len(n) - 1) %% block - street / 2) < 4,
                                             abs((seq_len(n) - 1) %% block - street / 2) < 4, `|`),
                    canopy)
canopy[buildings > 0] <- 0

dsm_city <- dtm + pmax(buildings, canopy)
green_city <- (canopy > 0.5 | park) * 1

# Open terrain variant: no buildings, sparse trees
canopy_open <- add_trees(4000, matrix(TRUE, n, n), matrix(0, n, n))
dsm_open <- dtm + canopy_open
green_open <- (canopy_open > 0.5 | smooth_field(n, k = 4, wl = c(200, 600)) > 0.2) * 1

to_rast <- function(m) {
  r <- rast(nrows = n, ncols = n, xmin = xmin, xmax = xmin + n, ymin = ymax - n, ymax = ymax, crs = crs)
  values(r) <- as.vector(t(m))
  r
}
writeRaster(to_rast(dsm_city), file.path(out_dir, "city_dsm.tif"), overwrite = TRUE)
writeRaster(to_rast(dtm), file.path(out_dir, "city_dtm.tif"), overwrite = TRUE)
writeRaster(to_rast(green_city), file.path(out_dir, "city_greenspace.tif"), overwrite = TRUE, datatype = "INT1U")
writeRaster(to_rast(dsm_open), file.path(out_dir, "open_dsm.tif"), overwrite = TRUE)
writeRaster(to_rast(dtm), file.path(out_dir, "open_dtm.tif"), overwrite = TRUE)
writeRaster(to_rast(green_open), file.path(out_dir, "open_greenspace.tif"), overwrite = TRUE, datatype = "INT1U")

# Observers: street centre lines every 10 m (inner part, away from the edges)
centre <- (seq(0, n - block, by = block) + street / 2)
centre <- centre[centre > 300 & centre < n - 300]
along <- seq(305, n - 305, by = 10)
pts <- rbind(
  expand.grid(row = centre, col = along),
  expand.grid(row = along, col = centre)
)
pts <- unique(pts)
obs_city <- st_as_sf(data.frame(x = xmin + pts$col - 0.5, y = ymax - pts$row + 0.5),
                     coords = c("x", "y"), crs = crs)
st_write(obs_city, file.path(out_dir, "observers_city.gpkg"), delete_dsn = TRUE, quiet = TRUE)

grid <- expand.grid(row = seq(310, n - 310, by = 20), col = seq(310, n - 310, by = 20))
obs_open <- st_as_sf(data.frame(x = xmin + grid$col - 0.5, y = ymax - grid$row + 0.5),
                     coords = c("x", "y"), crs = crs)
st_write(obs_open, file.path(out_dir, "observers_open.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# GAVI / lacunarity rasters
for (m in c(500, 1000, 2000)) {
  f1 <- smooth_field(m, k = 8, wl = c(20, 120))
  binary <- (f1 + matrix(rnorm(m * m, 0, 0.15), m, m) > 0.1) * 1
  cont <- pmin(pmax(0.5 + smooth_field(m, k = 8, wl = c(30, 200)) + matrix(rnorm(m * m, 0, 0.05), m, m), 0), 1)
  r <- c(rast(binary), rast(cont))
  names(r) <- c("greenspace", "ndvi")
  ext(r) <- c(0, m, 0, m)
  writeRaster(r, file.path(out_dir, sprintf("gavi_%d.tif", m)), overwrite = TRUE)
}

cat("Benchmark data written to", normalizePath(out_dir), "\n")
cat("City observers:", nrow(obs_city), " open-terrain observers:", nrow(obs_open), "\n")
cat(sprintf("Green share city: %.2f, open: %.2f\n", mean(green_city), mean(green_open)))
