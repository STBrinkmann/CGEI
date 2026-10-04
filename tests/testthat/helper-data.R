# Helpers to load the shipped test data and to call the C++ functions directly.

testdata_path <- function(...) testthat::test_path("testdata", ...)

# The test scene (see testdata/make_testdata.R).
load_scene <- function() {
  read <- function(f) {
    r <- terra::rast(testdata_path(f))
    terra::crs(r) <- "EPSG:25832"
    r
  }
  obs <- utils::read.csv(testdata_path("observers.csv"))
  list(
    dsm = read("scene_dsm.asc"),
    dtm = read("scene_dtm.asc"),
    greenspace = read("scene_greenspace.asc"),
    greenspace_fine = read("scene_greenspace_fine.asc"),
    observers = sf::st_as_sf(obs, coords = c("x", "y"), crs = 25832, remove = FALSE),
    expected = utils::read.csv(testdata_path("expected_scene.csv"))
  )
}

# SpatRaster from a [row, col] matrix (row 1 = north).
mk_rast <- function(mat, res = 1, xmin = 1000, ymax = 5000, crs = "EPSG:25832") {
  r <- terra::rast(nrows = nrow(mat), ncols = ncol(mat),
                   xmin = xmin, xmax = xmin + ncol(mat) * res,
                   ymin = ymax - nrow(mat) * res, ymax = ymax, crs = crs)
  terra::values(r) <- as.vector(t(mat))
  r
}

# Point observers (sf) at the centres of the given cells (1-based rows / cols).
mk_observers <- function(r, rows, cols) {
  xy <- terra::xyFromCell(r, terra::cellFromRowCol(r, rows, cols))
  sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"),
               crs = sf::st_crs(terra::crs(r)))
}

# Direct calls of the C++ functions. rows / cols are R's 1-based numbers.
cpp_vvi <- function(dsm, rows, cols, h0, radius, cores = 1L, early_stop = TRUE) {
  CGEI:::VVI_cpp(CGEI:::raster_geometry(dsm), terra::values(dsm, mat = FALSE),
                 as.integer(cols), as.integer(rows), h0, radius,
                 ncores = as.integer(cores), early_stop = early_stop)
}
cpp_vgvi <- function(dsm, green, rows, cols, h0, radius, fun = 3L, m = 1, b = 6,
                     cores = 1L, early_stop = TRUE) {
  CGEI:::VGVI_cpp(CGEI:::raster_geometry(dsm), terra::values(dsm, mat = FALSE),
                  CGEI:::raster_geometry(green), as.numeric(terra::values(green, mat = FALSE)),
                  as.integer(cols), as.integer(rows), h0, radius, as.integer(fun), m, b,
                  ncores = as.integer(cores), early_stop = early_stop)
}
cpp_rings <- function(dsm, green, rows, cols, h0, radius, cores = 1L, early_stop = TRUE) {
  out <- CGEI:::VGVI_rings_cpp(CGEI:::raster_geometry(dsm), terra::values(dsm, mat = FALSE),
                               CGEI:::raster_geometry(green),
                               as.numeric(terra::values(green, mat = FALSE)),
                               as.integer(cols), as.integer(rows), h0, radius,
                               ncores = as.integer(cores), early_stop = early_stop)
  lapply(out, as.data.frame)
}

# Random "urban" DSM: low noise plus random blocks of buildings / trees.
random_dsm <- function(nr, nc, n_blocks = 25, na_frac = 0, seed = 1) {
  set.seed(seed)
  m <- matrix(stats::runif(nr * nc, 0, 0.5), nr, nc)
  for (k in seq_len(n_blocks)) {
    r0 <- sample.int(nr, 1)
    c0 <- sample.int(nc, 1)
    h <- stats::runif(1, 2, 20)
    sz <- sample(1:4, 2, replace = TRUE)
    m[r0:min(nr, r0 + sz[1]), c0:min(nc, c0 + sz[2])] <- h
  }
  if (na_frac > 0) m[sample.int(nr * nc, round(na_frac * nr * nc))] <- NA
  m
}

has_openmp <- function() isTRUE(CGEI:::cgei_openmp_info()$openmp)
