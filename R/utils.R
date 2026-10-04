# Internal helpers shared by the R wrappers of the C++ functions.

# Raster geometry as expected by the C++ code (struct RasterInfo):
# c(xmin, xmax, ymin, ymax, nrow, ncol).
raster_geometry <- function(x) {
  e <- as.vector(terra::ext(x))
  c(e[["xmin"]], e[["xmax"]], e[["ymin"]], e[["ymax"]], terra::nrow(x), terra::ncol(x))
}

# Validate the `cores` argument and return it as integer. Warns (once per
# session) if more than one core is requested but CGEI was built without OpenMP.
check_cores <- function(cores) {
  checkmate::assert_count(cores, positive = TRUE, .var.name = "cores")
  cores <- as.integer(cores)
  if (cores > 1L && !isTRUE(cgei_openmp_info()$openmp) &&
      !isTRUE(getOption("CGEI.openmp_warned"))) {
    options(CGEI.openmp_warned = TRUE)
    warning("CGEI was built without OpenMP support; computations run on a single core.",
            call. = FALSE)
  }
  cores
}
