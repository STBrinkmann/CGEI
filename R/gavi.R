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
#' @importFrom terra classify
#' @keywords internal
reclassify_jenks <- function(raster_layer, n_classes = 9) {
  # Extract values from the raster layer (at most 50,000 of the non-NA cells)
  n_valid <- terra::global(raster_layer, "notNA")[1, 1]
  values <- unlist(terra::spatSample(raster_layer, min(50000, n_valid), na.rm = TRUE), use.names = FALSE)
  style <- ifelse(terra::ncell(raster_layer) > 5000, "fisher", "jenks")

  # Compute Jenks natural breaks
  breaks <- jenks_breaks(values, n_classes, style = style)
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf

  # Reclassify the raster layer based on the breaks
  # Create a matrix for reclassification, with the lower limit, upper limit, and new class value
  # (one row per class; there are fewer than n_classes classes if the layer
  # has fewer distinct values)
  rcl_mat <- cbind(utils::head(breaks, -1), utils::tail(breaks, -1), seq_len(length(breaks) - 1))
  reclassified_layer <- terra::classify(raster_layer, rcl_mat, include.lowest = TRUE)

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
  
  # Convert matrix to raster
  lac_mean_rast <- terra::rast(x)
  lac_mean_rast[] <- lac_mean_mat
  lac_mean_rast <- lac_mean_rast %>% 
    terra::crop(x, mask = TRUE)
  
  # Apply jenks on each layer to reclasify from 1-9
  if(progress) message("Reclassifying layers")
  for (i in 1:(terra::nlyr(lac_mean_rast))) {
    lac_mean_rast[[i]] <- reclassify_jenks(lac_mean_rast[[i]])
  }

  # Combine both rasters into one using mean
  gavi <- sum(lac_mean_rast) / terra::nlyr(lac_mean_rast)
  if(terra::nlyr(lac_mean_rast) > 1) {
    gavi <-  reclassify_jenks(gavi)
  }

  return(gavi)
}