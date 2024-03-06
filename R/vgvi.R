#' @title Viewshed Greenness Visibility Index (VGVI)
#'
#' @description The VGVI expresses the proportion of visible greenness to the total visible area based on a viewshed (\code{\link[CGEI]{viewshed_list}}).
#' The estimated VGVI values range between 0 and 1, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
#' A distance decay function is applied, to account for the reducing visual prominence of an object in space with increasing distance from the observer.
#'
#' @param observer An `sf` object containing the observer locations. Each observer location should be a point geometry with a defined coordinate reference system (CRS).
#' @param dsm_rast A `SpatRaster` object representing the Digital Surface Model (DSM) of the area, indicating the elevation of surface objects including vegetation, buildings, and other features.
#' @param dtm_rast A `SpatRaster` object representing the Digital Terrain Model (DTM) of the area, indicating the elevation of the ground surface without any objects.
#' @param greenspace_rast A `SpatRaster` object representing the greenspace within the area. Greenspaces are indicated by a value of 1, and non-greenspaces by a value of 0.
#' @param max_distance The maximum distance (in meters) from the observer location within which the greenspace visibility is calculated. Default is 200 meters.
#' @param observer_height The height (in meters) of the observer above the ground level. Default is 1.7 meters.
#' @param spacing optional; numeric > 0; Only if \code{observer} is a LINESTRING or POLYGON. Points on the line will be generated. The \code{spacing} parameter sets the distance in between the points on the line/grid. Defaults to the resolution of the \code{dsm_rast}.
#' @param m Parameter for the decay function, applicable if a decay function is selected. Default is 1.
#' @param b Parameter for the decay function, applicable if a decay function is selected. Default is 6.
#' @param mode A character string specifying the type of decay function to apply to the visibility weights. Options are "none" (no decay), "exponential", or "logit". Default is "none".
#' @param cores The number of cores to use for parallel processing. This parameter is relevant only if the function is set to run in parallel. Default is 1.
#' @param progress logical; Show progress bar and computation time?
#'
##' @details 
#' observer needs to be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON. If observer is a LINESTRING or MULTILINESTRING, 
#' points will be generated along the line(s) every "spacing" meters. If observer is a POLYGON or MULTIPOLYGON, a grid with resolution = "resolution" 
#' will be generated, and VGVI will be computed for every point.
#' The CRS (\code{\link[sf]{st_crs}}) needs to have a metric unit!
#' 
#' The type of function, used for calculating the distance decay weights, can be defined with the \code{mode} parameter.
#' The argument 'logit' uses the logistic function, d = 1 / (1 + e^(b * (x - m))) and 'exponential' the exponential function d = 1 / (1 + (b * x^m)).
#' The decay function can be visualized using the \code{\link[CGEI]{visualizeWeights}} function.
#'
#' @return sf_object containing the weighted VGVI values as POINT features, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' 
#' # Simulate observer locations as sf points
#' observer <- st_as_sf(data.frame(lon = c(5, 5.1), lat = c(52, 52.1)),
#'                       coords = c("lon", "lat"), crs = 4326)
#' 
#' # Transform to UTM zone 31N for metric units
#' observer <- st_transform(observer, 32631)
#' 
#' # Create synthetic rasters for DSM, DTM, and greenspace around observers
#' # This is a simplified example assuming flat terrain and random greenspaces
#' bbox_observer <- st_bbox(observer)
#' x_range <- bbox_observer[c("xmin", "xmax")]
#' y_range <- bbox_observer[c("ymin", "ymax")]
#' 
#' # Create a raster with 1km x 1km around the observer
#' dsm_rast <- rast(res=100, 
#'                  xmin=min(x_range)-1000, xmax=max(x_range) + 1000, 
#'                  ymin=min(y_range)-1000, ymax=max(y_range) + 1000, 
#'                  crs=crs(observer))
#' dsm_rast[] <- runif(ncell(dsm_rast), 0, 1.5) # Assign random heights
#' 
#' dtm_rast <- rast(dsm_rast, vals=0) # Flat terrain
#' 
#' greenspace_rast <- rast(dsm_rast)
#' values(greenspace_rast) <- sample(0:1, ncell(dsm_rast), replace=TRUE)
#' 
#' # Calculate VGVI
#' vgvi_results <- vgvi(observer, dsm_rast, dtm_rast, greenspace_rast)
#' print(vgvi_results)
#' }
#'
#' @export
#' @importFrom methods is
#' @importFrom sf st_crs st_as_sf st_transform st_geometry_type st_union st_cast st_line_sample st_set_geometry st_bbox st_buffer st_coordinates st_as_sfc
#' @importFrom dplyr rename mutate relocate everything
#' @importFrom terra crs rast res crop mask vect xyFromCell extract cellFromXY colFromX rowFromY writeRaster
#' @importFrom raster raster
#' @importFrom checkmate assert
#' @importFrom utils txtProgressBar setTxtProgressBar
vgvi <- function(observer, dsm_rast, dtm_rast, greenspace_rast,
                 max_distance = 200, observer_height = 1.7, 
                 spacing = NULL,
                 m = 1, b = 6, mode = c("none", "exponential", "logit"),
                 cores = 1L, progress = FALSE) {
  #### 1. Check input ####
  # Check observer
  valid_sf_types <- c("POINT", "MULTIPOINT", "LINESTRING", "MULTILINESTRING", "POLYGON", "MULTIPOLYGON")
  checkmate::assert(is(observer, "sf"), "observer must be an sf object")
  checkmate::assert(!sf::st_is_longlat(sf::st_crs(observer)), "observer must have a metric CRS")
  checkmate::assert(sf::st_geometry_type(observer, by_geometry = FALSE) %in% valid_sf_types, "observer must be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON")
  
  # Check rasters
  checkmate::assert(inherits(dsm_rast, "SpatRaster"), "dsm_rast must be a SpatRast object")
  checkmate::assert(inherits(dtm_rast, "SpatRaster"), "dtm_rast must be a SpatRast object")
  checkmate::assert(inherits(greenspace_rast, "SpatRaster"), "greenspace_rast must be a SpatRast object")
  
  # Check that the rasters have the same resolution
  checkmate::assert(length(unique(terra::res(dsm_rast))) == 1, "dsm_rast: x and y resolution must be equal")
  checkmate::assert(all(terra::res(dsm_rast) == terra::res(dtm_rast)), "dsm_rast and dtm_rast must have the same resolution")
  
  # Check crs
  checkmate::assert(identical(sf::st_crs(observer)$proj4string, terra::crs(dtm_rast, proj = TRUE)), "dsm_rast and dtm_rast must have the same crs")
  checkmate::assert(identical(sf::st_crs(observer)$proj4string, terra::crs(greenspace_rast, proj = TRUE)), "dsm_rast and greenspace_rast must have the same crs")
  checkmate::assert(identical(sf::st_crs(observer)$proj4string, terra::crs(dsm_rast, proj = TRUE)), "dsm_rast and observer must have the same crs")
  
  # Check max_distance
  checkmate::assert(methods::is(max_distance, "numeric"), "max_distance must be a numeric")
  checkmate::assert(max_distance > 0, "max_distance must be greater than 0")
  max_distance <- round(max_distance, digits = 0)
  
  # Check spacing
  checkmate::assert(methods::is(spacing, "numeric") | is.null(spacing), "spacing must be a numeric or NULL")
  if(is.null(spacing)) {
    spacing <- terra::res(dsm_rast)[1]
  }
  
  # Check observer_height
  checkmate::assert(methods::is(observer_height, "numeric"), "observer_height must be a numeric")
  checkmate::assert(observer_height > 0, "observer_height must be greater than 0")
  
  # Check m and b
  checkmate::assert(methods::is(m, "numeric"), "m must be a numeric")
  checkmate::assert(methods::is(b, "numeric"), "b must be a numeric")
  
  # Check mode
  mode <- match.arg(mode, c("none", "exponential", "logit"))
  mode <- switch(mode,
                 logit = 1,
                 exponential = 2,
                 none = 3)
  
  
  #### 2. Convert observer to points
  if(progress) {
    message("Preprocessing:")
    pb = utils::txtProgressBar(min = 0, max = 5, initial = 0, style = 2)
  }
  observer <- sf_to_POINT(observer, spacing, dsm_rast)
  if (progress) utils::setTxtProgressBar(pb, 1)
  
  
  #### 3. Prepare data for viewshed analysis ####
  # Max AOI
  max_aoi <- observer %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_buffer(max_distance)
  
  # Crop DSM to max AOI
  dsm_rast <- terra::crop(dsm_rast, terra::vect(max_aoi))
  greenspace_rast <- terra::crop(greenspace_rast, terra::vect(max_aoi))
  
  dsm_vec <- terra::values(dsm_rast, mat = FALSE)
  greenspace_vec <- terra::values(greenspace_rast, mat = FALSE)
  
  dsm_cpp_rast <- dsm_rast %>% terra::rast() %>% raster::raster()
  greenspace_cpp_rast <- greenspace_rast %>% terra::rast() %>% raster::raster()
  
  # Coordinates of start point
  x0 <- sf::st_coordinates(observer)[,1]
  y0 <- sf::st_coordinates(observer)[,2]
  
  # Observer heights
  height_0_vec <- unlist(terra::extract(dtm_rast, cbind(x0, y0)), use.names = F) + observer_height
  
  if (progress) utils::setTxtProgressBar(pb, 2)
  
  
  #### 4. Remove points outside the DSM or DTM ####
  invalid_points <- unique(c(
    which(is.na(terra::extract(dsm_rast, cbind(x0, y0)))), # points outside the DSM
    which(is.na(height_0_vec)) # points outside the DTM
  ))
  
  # Remove invalid points
  if (length(invalid_points) > 0) {
    observer <- observer[-invalid_points, ]
    x0 <- x0[-invalid_points]
    y0 <- y0[-invalid_points]
    height_0_vec <- height_0_vec[-invalid_points]
  }
  
  if (progress) utils::setTxtProgressBar(pb, 3)
  
  
  #### 5. Last steps of PreProcessing ####
  # Prepare observer for output
  observer <- observer %>% 
    dplyr::mutate(VGVI = as.numeric(NA)) %>% 
    dplyr::select(VGVI, dplyr::everything())
  
  # convert x0/y0 to col/row
  c0 <- terra::colFromX(dsm_rast, x0)
  r0 <- terra::rowFromY(dsm_rast, y0)
  
  if (progress) utils::setTxtProgressBar(pb, 4)
  
  
  # 6. Final update of Pre-Processing ProgressBar
  if (progress) {
    utils::setTxtProgressBar(pb, 5)
    cat("\n")
  }
  if (length(invalid_points) == 1) {
    message("1 point has been removed, because it was outside of the DSM or DTM")
  } else if (length(invalid_points) > 1) {
    message(paste(length(invalid_points), "points have been removed, because they were outside of the DSM or DTM"))
  }
  invisible(gc())
  
  
  #### 7. Calculate viewsheds and VGVI ####
  if (progress) {
    message(paste0("Computing VGVI for ", nrow(observer), ifelse(nrow(observer)>1, " points:", " point:")))
    start_time <- Sys.time()
    
    on.exit({
      # Runtime statistics
      time_dif <- round(cores * ((as.numeric(difftime(Sys.time(), start_time, units = "s"))*1000) / nrow(observer)), 2)
      cat("\n")
      
      time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "m")))
      if(time_total < 1){
        time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "s")))
        
        if(time_total < 1){
          time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "s")))*1000
          message(paste("Total runtime:", time_total, " milliseconds"))
        } else {
          message(paste("Total runtime:", time_total, " seconds"))
        }
      } else {
        message(paste("Total runtime:", time_total, " minutes"))
      }
      
      message(paste("Average time for a single point:", time_dif, "milliseconds"))
    })
  }
      
  
  vgvi_values <- VGVI_cpp(dsm = dsm_cpp_rast, dsm_values = dsm_vec,
                          greenspace = greenspace_cpp_rast, greenspace_values = greenspace_vec,
                          x0 = c0, y0 = r0, radius = max_distance, h0 = height_0_vec,
                          fun = mode, m = m, b = b, ncores = cores, display_progress = progress)
  
  valid_values <- unlist(lapply(vgvi_values, is.numeric), use.names = FALSE)
  observer[["VGVI"]][valid_values] <- vgvi_values[valid_values]
  
  return(observer)
}