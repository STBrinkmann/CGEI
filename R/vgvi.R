#' Calculate Weights Based on Distance and Decay Mode
#'
#' This internal function calculates weights for each cell in a distance raster based on the specified decay mode. The function supports logistic and exponential decay modes, applying a decay function to the normalized distance to compute weights. This function is typically used internally within spatial analysis functions to model the diminishing influence of features (e.g., greenspaces) with increasing distance from a point of interest.
#'
#' @param distance_raster A `SpatRaster` object representing distances from a specific point or points. Each cell's value represents the distance to the point(s) of interest.
#' @param mode A character string specifying the type of decay function to apply. Supported modes are "logistic" and "exponential". If any other value is provided, the function returns a weight of 1 for all distances, implying no decay.
#' @param m A numeric value used as a parameter in the decay functions. For the logistic mode, `m` adjusts the midpoint of the sigmoid curve where the weight is 0.5. For the exponential mode, `m` influences the rate at which weights decrease with distance.
#' @param b A numeric value used as a parameter in the decay functions to control the steepness (logistic) or the scale (exponential) of the decay.
#'
#' @return A `SpatRaster` object where each cell's value represents the calculated weight based on the specified decay mode and the distance to the point(s) of interest.
#'
#' @keywords internal
calculate_weights <- function(distance_raster, mode, m, b) {
  # Normalize the distance raster
  norm_distance <- distance_raster / max(distance_raster[], na.rm = TRUE)
  
  # Select the decay function based on the input
  if (mode == 'logistic') {
    return(weight_raster <- 1 / (1 + exp(b * (norm_distance - m))))
  } else if (mode == 'exponential') {
    return(weight_raster <- 1 / (1 + (b * (norm_distance ^ m))))
  } else {
    return(1)
  }
}


#' Viewshed Greenness Visibility Index (VGVI)
#'
#' @description The VGVI expresses the proportion of visible greenness to the total visible area based on a \code{\link[terra]{viewshed}}.
#' The estimated VGVI values range between 0 and 1, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
#' NA values are returned for observer locations where the observer height is greater than the maximum elevation in the DSM raster.
#' A distance decay function can be applied, to account for the reducing visual prominence of an object in space with increasing distance from the observer.
#'
#' @param observers An `sf` object containing the observer locations. Each observer location should be a point geometry with a defined coordinate reference system (CRS).
#' @param dsm_rast A `SpatRaster` object representing the Digital Surface Model (DSM) of the area, indicating the elevation of surface objects including vegetation, buildings, and other features.
#' @param dtm_rast A `SpatRaster` object representing the Digital Terrain Model (DTM) of the area, indicating the elevation of the ground surface without any objects.
#' @param greenspace_rast A `SpatRaster` object representing the greenspace within the area. Greenspaces are indicated by a value of 1, and non-greenspaces by a value of 0.
#' @param max_distance The maximum distance (in meters) from the observer location within which the greenspace visibility is calculated. Default is 200 meters.
#' @param observer_height The height (in meters) of the observer above the ground level. Default is 1.7 meters.
#' @param spacing optional; numeric > 0; Only if \code{observer} is a LINESTRING. Points on the line will be generated. The \code{spacing} parameter sets the distance in between the points on the line/grid. Defaults to the resolution of the \code{dsm_rast}.
#' @param m Parameter for the decay function, applicable if a decay function is selected. Default is 1.
#' @param b Parameter for the decay function, applicable if a decay function is selected. Default is 6.
#' @param mode A character string specifying the type of decay function to apply to the visibility weights. Options are "none" (no decay), "exponential", or "logit". Default is "none".
#' @param cores The number of cores to use for parallel processing. This parameter is relevant only if the function is set to run in parallel. Default is 1.
#'
#' @return A numeric vector containing the VGVI values for each observer location.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' 
#' # Simulate observer locations as sf points
#' observers <- st_as_sf(data.frame(lon = c(5, 5.1), lat = c(52, 52.1)),
#'                       coords = c("lon", "lat"), crs = 4326)
#' 
#' # Transform to UTM zone 31N for metric units
#' observers <- st_transform(observers, 32631)
#' 
#' # Create synthetic rasters for DSM, DTM, and greenspace around observers
#' # This is a simplified example assuming flat terrain and random greenspaces
#' bbox_observers <- st_bbox(observers)
#' x_range <- bbox_observers[c("xmin", "xmax")]
#' y_range <- bbox_observers[c("ymin", "ymax")]
#' 
#' # Create a raster with 1km x 1km around the observers
#' dsm_rast <- rast(nrows=100, ncols=100,
#'                  xmin=min(x_range)-1000, xmax=max(x_range) + 1000, 
#'                  ymin=min(y_range)-1000, ymax=max(y_range) + 1000, 
#'                  crs=crs(observers))
#' dsm_rast[] <- runif(ncell(dsm_rast), 0, 1.5) # Assign random heights
#' 
#' dtm_rast <- rast(dsm_rast, vals=0) # Flat terrain
#' 
#' greenspace_rast <- rast(dsm_rast)
#' values(greenspace_rast) <- sample(0:1, ncell(dsm_rast), replace=TRUE)
#' 
#' # Calculate VGVI
#' vgvi_results <- vgvi(observers, dsm_rast, dtm_rast, greenspace_rast)
#' print(vgvi_results)
#' }
#'
#' @export
#' @importFrom methods is
#' @importFrom sf st_as_sf st_transform st_geometry_type st_union st_cast st_line_sample st_set_geometry
#' @importFrom terra rast res crop extract writeRaster buffer viewshed distance
#' @importFrom checkmate assert
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
vgvi <- function(observers, dsm_rast, dtm_rast, greenspace_rast,
                 max_distance = 200, observer_height = 1.7, 
                 spacing = NULL,
                 m = 1, b = 6, mode = c("none", "exponential", "logit"),
                 cores = 1L) {
  optProgress <- terra::terraOptions(print = FALSE)$progress
  terra::terraOptions(progress=0)
  
  # Check observers
  checkmate::assert(is(observers, "sf"), "observers must be an sf object")
  
  # Check rasters
  checkmate::assert(inherits(dsm_rast, "SpatRaster"), "dsm_rast must be a SpatRast object")
  checkmate::assert(inherits(dtm_rast, "SpatRaster"), "dtm_rast must be a SpatRast object")
  checkmate::assert(inherits(greenspace_rast, "SpatRaster"), "greenspace_rast must be a SpatRast object")
  
  # Check that the rasters have the same resolution
  checkmate::assert(all(terra::res(dsm_rast) == terra::res(dtm_rast)), "dsm_rast and dtm_rast must have the same resolution")
  checkmate::assert(all(terra::res(dsm_rast) == terra::res(greenspace_rast)), "dsm_rast and greenspace_rast must have the same resolution")
  
  # Check other parameters
  checkmate::assert(methods::is(max_distance, "numeric"), "max_distance must be a numeric")
  checkmate::assert(methods::is(observer_height, "numeric"), "observer_height must be a numeric")
  checkmate::assert(methods::is(m, "numeric"), "m must be a numeric")
  checkmate::assert(methods::is(b, "numeric"), "b must be a numeric")
  checkmate::assert(methods::is(mode, "character"), "mode must be a character")
  mode <- match.arg(mode, c("none", "exponential", "logit"))
  
  # Check crs
  checkmate::assert(identical(sf::st_crs(observers)$proj4string, terra::crs(dtm_rast, proj = TRUE)), "dsm_rast and dtm_rast must have the same crs")
  checkmate::assert(identical(sf::st_crs(observers)$proj4string, terra::crs(greenspace_rast, proj = TRUE)), "dsm_rast and greenspace_rast must have the same crs")
  checkmate::assert(identical(sf::st_crs(observers)$proj4string, terra::crs(dsm_rast, proj = TRUE)), "dsm_rast and observer must have the same crs")
  
  # Check spacing
  if(!is.null(spacing)) {
    checkmate::assert(methods::is(spacing, "numeric"), "spacing must be a numeric")
  } else {
    spacing <- terra::res(dsm_rast)[1]
  }
  
  # Convert observer to POINT if it is LINESTRING/MULTILINESTRING
  if(sf::st_geometry_type(observers, by_geometry = FALSE) %in% c("LINESTRING", "MULTILINESTRING")) {
    geom_name <- 
    observers <- observers %>%
      sf::st_union() %>% 
      sf::st_cast("LINESTRING") %>%
      sf::st_line_sample(density = 1/2)
    
    observers <- observers[!sf::st_is_empty(observers)] %>% 
      sf::st_cast("POINT") %>% 
      sf::st_as_sf() %>% 
      sf::st_set_geometry("geom")
  }
  
  # Calculate observer height above ground
  observers_height_alg <- observer_height + 
    as.numeric(unlist(terra::extract(dtm_rast, observers, ID = FALSE))) - 
    as.numeric(unlist(terra::extract(dsm_rast, observers, ID = FALSE)))
  
  # Save the raster in temporary files to be used in parallel
  if(cores > 1) {
    dsm_rast_p <- tempfile(fileext = ".tif")
    dtm_rast_p <- tempfile(fileext = ".tif")
    greenspace_rast_p <- tempfile(fileext = ".tif")
    terra::writeRaster(dsm_rast, file = dsm_rast_p)
    terra::writeRaster(dtm_rast, file = dtm_rast_p)
    terra::writeRaster(greenspace_rast, file = greenspace_rast_p)
  }
  
  # Set up the future plan for parallel backend
  future::plan(future::multisession, workers = cores)
  
  # Process each observer location in parallel
  vgv_indices <- future.apply::future_lapply(1:nrow(observers), 
                                             future.packages = c("terra", "sf"),
                                             future.seed=TRUE,
                                             function(i) {
                                               if(cores > 1) {
                                                 dsm_rast <- terra::rast(dsm_rast_p)
                                                 dtm_rast <- terra::rast(dtm_rast_p)
                                                 greenspace_rast <- terra::rast(greenspace_rast_p)
                                               }
                                               
                                               # Current observer
                                               observer <- observers[i, , drop = FALSE]
                                               observer_height_alg <- observers_height_alg[i]
                                               
                                               if(observer_height_alg < 0) {
                                                 warning("Observer height is below ground level. Returning NA.")
                                                 return(NA)
                                               }
                                               # Get a buffer mask for the observer with the max distance
                                               observer_buffer <- terra::buffer(terra::vect(observer), max_distance)
                                               
                                               # Mask the rasters (DSM and GreenSpace)
                                               dsm_rast_aoi <- terra::crop(dsm_rast, observer_buffer, mask = TRUE)
                                               greenspace_rast_aoi <- terra::crop(greenspace_rast, observer_buffer, mask = F)
                                               
                                               # Estimate viewshed
                                               viewshed <- terra::viewshed(x = dsm_rast_aoi, 
                                                                           loc = sf::st_coordinates(observer), 
                                                                           observer = observer_height_alg)
                                               
                                               # Calculate the distance from the observer to each cell
                                               distance <- terra::distance(x = viewshed, y = terra::vect(observer))
                                               
                                               # Calculate the visible green cells (viewshed == TRUE & greenspace_rast_aoi == 1)
                                               visible_green <- viewshed & greenspace_rast_aoi == 1
                                               
                                               # Apply the selected decay function to calculate weights
                                               weights <- calculate_weights(distance, mode, m, b)
                                               
                                               # Create weighted rasters
                                               weighted_green <- visible_green * weights
                                               weighted_viewshed <- viewshed * weights
                                               
                                               # Calculate VGVI with weights for the current observer
                                               return(sum(weighted_green[], na.rm = TRUE) / sum(weighted_viewshed[], na.rm = TRUE))
                                             })
  
  # Clean up temporary files
  if(cores > 1) {
    file.remove(dsm_rast_p)
    file.remove(dtm_rast_p)
    file.remove(greenspace_rast_p)
  }
  terra::terraOptions(progress=optProgress, print=FALSE)
  
  observers <- observers %>% 
    dplyr::mutate(VGVI = unlist(vgv_indices)) %>%
    dplyr::relocate(VGVI)
  
  return(observers)
}