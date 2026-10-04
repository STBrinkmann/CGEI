#' Natural breaks
#'
#' Breaks of the natural breaks classification of `values` into `n_classes`
#' classes, identical to `classInt::classIntervals(values, n_classes, style)$brks`
#' for `style = "jenks"` and the same optimal partition for `style = "fisher"`
#' (see `jenks_breaks_cpp()`), but orders of magnitude faster.
#' With fewer distinct values than classes, every distinct value becomes a class
#' (classInt's "unique" breaks).
#'
#' @noRd
#' @keywords internal
jenks_breaks <- function(values, n_classes = 9, style = c("fisher", "jenks")) {
  style <- match.arg(style)
  values <- values[is.finite(values)]
  u <- sort(unique(values))
  if (length(u) == 0) stop("No finite values to classify.")
  if (length(u) == 1) return(c(u, u))
  n <- min(as.integer(n_classes), length(u))
  if (n == length(u)) {
    du <- diff(u)
    return(c(u[1] - mean(du) / 2, u[-length(u)] + du / 2, u[length(u)] + mean(du) / 2))
  }
  jenks_breaks_cpp(values, n, style)
}

#' Classify values using Jenks natural breaks
#'
#' Classifies `values` into (at most) `n_classes` natural breaks classes 1, 2, ...
#' The breaks are computed from all non-NA values or, if there are more than
#' `max_sample`, from a random sample of `max_sample` of them. The classes are
#' closed on the right (the lowest one also on the left), as with
#' `terra::classify(..., include.lowest = TRUE)`; NA stays NA.
#'
#' @noRd
#' @param values numeric vector.
#' @param n_classes The number of classes; default = 9.
#' @param style "fisher" or "jenks", see `jenks_breaks()`.
#' @param max_sample Maximum number of values used to compute the breaks.
#'
#' @return integer vector of classes.
#' @keywords internal
classify_jenks <- function(values, n_classes = 9, style = "fisher", max_sample = 50000) {
  valid <- values[!is.na(values)]
  if (length(valid) > max_sample) valid <- valid[sample.int(length(valid), max_sample)]

  breaks <- jenks_breaks(valid, n_classes, style = style)
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  findInterval(values, breaks, left.open = TRUE, rightmost.closed = TRUE)
}

#' Reclassify Raster Layer Using Jenks Natural Breaks
#'
#' This function reclassifies a given raster layer into specified number of classes
#' based on the Jenks natural breaks classification method. It is designed to handle
#' large datasets by sampling up to 50,000 non-NA values from the raster layer for
#' computing the breaks. This method is particularly useful for categorizing continuous
#' data into natural clusters.
#'
#' @noRd
#' @param raster_layer A `SpatRaster` object to be reclassified.
#' @param n_classes The number of classes to divide the raster layer into; default=9.
#'
#' @return A reclassified `SpatRaster` object, with values categorized into the
#' specified number of classes based on the Jenks natural breaks.
#'
#' @keywords internal
reclassify_jenks <- function(raster_layer, n_classes = 9) {
  style <- if (terra::ncell(raster_layer) > 5000) "fisher" else "jenks"
  reclassified_layer <- terra::rast(raster_layer)
  terra::values(reclassified_layer) <- classify_jenks(terra::values(raster_layer, mat = FALSE),
                                                      n_classes, style = style)
  return(reclassified_layer)
}


#' Greenspace Availability Index (GAVI)
#'
#' This function computes the Greenspace Availability Index (GAVI) from a lacunarity dataset. 
#'
#' @param x A `SpatRaster` object.
#' @param lac A data frame containing lacunarity information with columns (see (\code{\link[CGEI]{lacunarity}}).
#' @param na.rm A logical indicating whether NA values should be removed.
#' @param cores The number of cores to use for parallel processing. Default is 1.
#' @param progress logical; Show progress bar?
#' 
#' @details
#' For every layer, the focal means of all box sizes in \code{lac} are weighted by their lacunarity
#' and averaged. These values are classified into 9 classes using Jenks natural breaks ("fisher"
#' style for rasters with more than 5,000 cells, "jenks" otherwise). The GAVI is the mean of the
#' classes of all layers, again classified into 9 classes if \code{x} has more than one layer.
#' Cells that are NA in a layer are NA in that layer's classes (and in the GAVI).
#' The natural breaks are computed from all non-NA cells or, if there are more than 50,000 of them,
#' from a random sample of 50,000 cells; use \code{set.seed()} for reproducible results.
#' 
#' @return A `SpatRaster` object representing the GAVI.
#'
#' @examples
#' library(CGEI)
#' library(terra)
#' 
#' mat_sample <- matrix(data = c(
#'   1,1,0,1,1,1,0,1,0,1,1,0,
#'   0,0,0,0,0,1,0,0,0,1,1,1,
#'   0,1,0,1,1,1,1,1,0,1,1,0,
#'   1,0,1,1,1,0,0,0,0,0,0,0,
#'   1,1,0,1,0,1,0,0,1,1,0,0,
#'   0,1,0,1,1,0,0,1,0,0,1,0,
#'   0,0,0,0,0,1,1,1,1,1,1,1,
#'   0,1,1,0,0,0,1,1,1,1,0,0,
#'   0,1,1,1,0,1,1,0,1,0,0,1,
#'   0,1,0,0,0,0,0,0,0,1,1,1,
#'   0,1,0,1,1,1,0,1,1,0,1,0,
#'   0,1,0,0,0,1,0,1,1,1,0,1
#' ), nrow = 12, ncol = 12, byrow = TRUE)
#' 
#' x <- rast(mat_sample)
#' x2 <- rast(mat_sample*runif(144, 1, 2))
#' 
#' x <- c(x, x2)
#' lac <- lacunarity(x)
#' gavi(x, lac)
#'
#' @importFrom terra values rast
#' @importFrom checkmate assert_class assert_set_equal assert_true
#' @export
gavi <- function(x, lac, na.rm = TRUE, cores = 1, progress = FALSE) {
  # Check input
  checkmate::assert_class(x, "SpatRaster")
  checkmate::assert_set_equal(names(lac), c("name", "i", "r", "ln(r)", "Lac", "ln(Lac)"))
  checkmate::assert_class(na.rm, "logical")
  checkmate::assert_true(length(unique(lac[["i"]])) == terra::nlyr(x))
  cores <- check_cores(cores)
  
  # Convert raster to matrix
  x_mat <- terra::values(x, mat = TRUE)
  storage.mode(x_mat) <- "double"
  x_rast <- raster_geometry(x)
  
  # Apply focal C++ function
  lac <- lac[,c("i", "r", "Lac")] %>% as.matrix()
  
  lac_mean_mat <- focal_sum(x = x_rast, x_mat = x_mat, lac = lac, na_rm = na.rm,
                            ncores = cores, display_progress = progress)
  # Cells that are NA in x are NA in the result (layer by layer)
  lac_mean_mat[is.na(x_mat)] <- NA
  rm(x_mat)
  
  # Apply jenks on each layer to reclasify from 1-9
  # (on the values in memory, the breaks from at most 50,000 sampled cells)
  if(progress) message("Reclassifying layers")
  style <- if (terra::ncell(x) > 5000) "fisher" else "jenks"
  for (i in seq_len(ncol(lac_mean_mat))) {
    lac_mean_mat[, i] <- classify_jenks(lac_mean_mat[, i], 9, style = style)
  }

  # Combine both layers into one using mean (NA if a layer is NA)
  gavi_vec <- rowSums(lac_mean_mat) / ncol(lac_mean_mat)
  if(ncol(lac_mean_mat) > 1) {
    gavi_vec <- classify_jenks(gavi_vec, 9, style = style)
  }

  gavi <- terra::rast(x, nlyrs = 1)
  terra::values(gavi) <- gavi_vec
  names(gavi) <- "sum"

  return(gavi)
}
