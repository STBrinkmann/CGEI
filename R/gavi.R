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
#' @importFrom classInt classIntervals
#' @keywords internal
reclassify_jenks <- function(raster_layer, n_classes = 9) {
  # Extract values from the raster layer
  values <- unlist(terra::spatSample(raster_layer, min(50000, terra::ncell(raster_layer)), na.rm = TRUE), use.names = FALSE)
  style <- ifelse(terra::ncell(raster_layer) > 5000, "fisher", "jenks")
  
  # Compute Jenks natural breaks
  breaks <- classInt::classIntervals(values, n_classes, style = style,
                                     warnLargeN = FALSE)$brks
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  
  # Reclassify the raster layer based on the breaks
  # Create a matrix for reclassification, with the lower limit, upper limit, and new class value
  rcl_mat <- matrix(c(head(breaks, -1), tail(breaks, -1), 1:n_classes), ncol = 3)
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
#' @importFrom raster raster
#' @importFrom checkmate assert_class assert_set_equal assert_true
#' @export
gavi <- function(x, lac, na.rm = TRUE, cores = 1, progress = FALSE) {
  # Check input
  checkmate::assert_class(x, "SpatRaster")
  checkmate::assert_set_equal(names(lac), c("name", "i", "r", "ln(r)", "Lac", "ln(Lac)"))
  checkmate::assert_class(na.rm, "logical")
  checkmate::assert_true(length(unique(lac[["i"]])) == terra::nlyr(x))
  
  # Convert raster to matrix
  x_mat <- terra::values(x, mat = TRUE)
  x_rast <- x %>% terra::rast() %>% raster::raster()
  
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