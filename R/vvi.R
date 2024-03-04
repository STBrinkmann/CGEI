#' Viewshed Visibility Index (VVI) from sf
#' 
#' @description The VVI expresses the proportion of visible area to the total area based on a viewshed.
#' The estimated VVI values range between 0 and 1, where 0 = no cells within the buffer is visible, and 1 = all of the cells are visible.
#'
#' @param observer An `sf` object containing the observer locations. Each observer location should be a point geometry with a defined coordinate reference system (CRS).
#' @param dsm_rast A `SpatRaster` object representing the Digital Surface Model (DSM) of the area, indicating the elevation of surface objects including vegetation, buildings, and other features.
#' @param dtm_rast A `SpatRaster` object representing the Digital Terrain Model (DTM) of the area, indicating the elevation of the ground surface without any objects.
#' @param max_distance The maximum distance (in meters) from the observer location within which the greenspace visibility is calculated. Default is 200 meters.
#' @param observer_height The height (in meters) of the observer above the ground level. Default is 1.7 meters.
#' @param spacing optional; numeric > 0; Only if \code{observer} is a LINESTRING or POLYGON. Points on the line will be generated. The \code{spacing} parameter sets the distance in between the points on the line/grid. Defaults to the resolution of the \code{dsm_rast}.
#' @param mode A character string specifying the type of decay function to apply to the visibility weights. Options are "VVI", or "viewshed". Default is "VVI".
#' @param cores The number of cores to use for parallel processing. This parameter is relevant only if the function is set to run in parallel. Default is 1.
#' @param progress logical; Show progress bar and computation time?
#'
##' @details 
#' observer needs to be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON. If observer is a LINESTRING or MULTILINESTRING, 
#' points will be generated along the line(s) every "resolution" meters. If observer is a POLYGON or MULTIPOLYGON, a grid with resolution = "resolution" 
#' will be generated, and VVI will be computed for every point.
#' The CRS (\code{\link[sf]{st_crs}}) needs to have a metric unit!
#' 
#' If \code{mode = "VVI"}, an sf_object containing the VVI values as POINT features, where 0 = no visible cells, and 1 = all of the cells are visible.
#' When mode is set to \code{viewshed}, the function outputs a SpatRaster with two layers derived from dsm_rast:
#' 1. \code{n_views}: Counts how many times each cell is visible across all viewsheds. This layer identifies cells with high visibility across the landscape.
#' 2. \code{view_per_viewshed}: Calculates the ratio of the number of viewsheds in which a cell is visible to the total number of potential viewsheds for that cell. This value, ranging between 0 and 1, indicates the proportion of potential observation points from which the cell is visible, providing a normalized measure of visibility.
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
#' # Calculate VGVI
#' vvi_results <- vvi(observer, dsm_rast, dtm_rast)
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
vvi <- function(observer, dsm_rast, dtm_rast,
                max_distance = 200, observer_height = 1.7, 
                spacing = NULL, mode = c("VVI", "viewshed"), 
                cores = 1L, progress = FALSE) {
  #### 1. Check input ####
  # Check observer
  checkmate::assert(is(observer, "sf"), "observer must be an sf object")
  checkmate::assert(!sf::st_is_longlat(sf::st_crs(observer)), "observer must have a metric CRS")
  
  # Check rasters
  checkmate::assert(inherits(dsm_rast, "SpatRaster"), "dsm_rast must be a SpatRast object")
  checkmate::assert(inherits(dtm_rast, "SpatRaster"), "dtm_rast must be a SpatRast object")
  
  # Check that the rasters have the same resolution
  checkmate::assert(length(unique(terra::res(dsm_rast))) == 1, "dsm_rast: x and y resolution must be equal")
  checkmate::assert(all(terra::res(dsm_rast) == terra::res(dtm_rast)), "dsm_rast and dtm_rast must have the same resolution")
  
  # Check crs
  checkmate::assert(identical(sf::st_crs(observer)$proj4string, terra::crs(dtm_rast, proj = TRUE)), "dsm_rast and dtm_rast must have the same crs")
  checkmate::assert(identical(sf::st_crs(observer)$proj4string, terra::crs(dsm_rast, proj = TRUE)), "dsm_rast and observer must have the same crs")
  
  # Check max_distance
  checkmate::assert(methods::is(max_distance, "numeric"), "max_distance must be a numeric")
  checkmate::assert(max_distance > 0, "max_distance must be greater than 0")
  max_distance <- round(max_distance, digits = 0)
  
  # Check observer_height
  checkmate::assert(methods::is(observer_height, "numeric"), "observer_height must be a numeric")
  checkmate::assert(observer_height > 0, "observer_height must be greater than 0")
  
  # Check spacing
  checkmate::assert(methods::is(spacing, "numeric") | is.null(spacing), "spacing must be a numeric or NULL")
  if(is.null(spacing)) {
    spacing <- terra::res(dsm_rast)[1]
  }

  # Check mode
  mode <- match.arg(mode, c("VVI", "viewshed"))
  
  #### 2. Convert observer to points
  if(progress) {
    message("Preprocessing:")
    pb = utils::txtProgressBar(min = 0, max = 5, initial = 0, style = 2)
  }
  observer <- sf_to_POINT(observer)
  
  if (progress) utils::setTxtProgressBar(pb, 1)
  
  
  #### 3. Prepare data for viewshed analysis ####
  # Max AOI
  max_aoi <- observer %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_buffer(max_distance)
  
  # Crop DSM to max AOI
  dsm_rast <- terra::crop(dsm_rast, terra::vect(max_aoi))
  dsm_vec <- terra::values(dsm_rast, mat = FALSE)
  dsm_cpp_rast <- dsm_rast %>% terra::rast() %>% raster::raster()
  
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
    message(paste0("Computing VVI for ", nrow(observer), ifelse(nrow(observer)>1, " points:", " point:")))
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

  vvi_list <- VVI_cpp(dsm = dsm_cpp_rast, dsm_values = dsm_vec,
                      x0 = c0, y0 = r0, radius = max_distance, h0 = height_0_vec,
                      ncores = cores, display_progress = progress)
  
  if(mode == "VVI") {
    visible_cells <- unlist(sapply(vvi_list, function(vvi) {
      length(vvi$visible_cells)
    }))
    viewshed <- unlist(sapply(vvi_list, function(vvi) {
      length(vvi$viewshed)
    }))
    vvi <- visible_cells/viewshed
    
    
    observer <- observer %>% 
      dplyr::mutate(VVI = vvi,
                    n_visible_cells = visible_cells,
                    n_total_cells = viewshed) %>% 
      dplyr::select(VVI, n_visible_cells, n_total_cells, dplyr::everything())
    return(observer)
  } else if(mode == "viewshed") {
    viewshed_cells <- unlist(sapply(vvi_list, function(vvi) {
      vvi$viewshed
    }))
    visible_cells <- unlist(sapply(vvi_list, function(vvi) {
      vvi$visible_cells
    }))
    
    viewshed_count <- tabulate(viewshed_cells, nbins = terra::ncell(dsm_rast))
    visible_count <- tabulate(visible_cells, nbins = terra::ncell(dsm_rast))
    
    output <- terra::rast(dsm_rast)
    output$n_views <- visible_count
    output[[1]][is.na(visible_count / viewshed_count)] <- NA
    output$view_per_viewshed <- visible_count / viewshed_count
    return(output)
  }
}

# Helper function that converts SF to POINT
sf_to_POINT <- function(x) {
  if (as.character(sf::st_geometry_type(x, by_geometry = FALSE)) %in% c("LINESTRING", "MULTILINESTRING")) {
    sf_column <- attr(x, "sf_column")
    x <- x %>%
      sf::st_union() %>% 
      sf::st_cast("LINESTRING") %>%
      sf::st_line_sample(density = 1/spacing)
    
    x <- x[!sf::st_is_empty(x)] %>% 
      sf::st_cast("POINT") %>% 
      sf::st_as_sf() %>% 
      sf::st_set_geometry(sf_column)
    rm(sf_column)
  } else if (as.character(sf::st_geometry_type(x, by_geometry = FALSE)) %in% c("POLYGON", "MULTIPOLYGON")) {
    sf_column <- attr(x, "sf_column")
    x_bbox <- sf::st_bbox(x)
    x <- terra::rast(xmin = x_bbox[1], xmax = x_bbox[3], 
                            ymin = x_bbox[2], ymax = x_bbox[4], 
                            crs = terra::crs(dsm_rast), resolution = spacing, vals = 0) %>% 
      terra::crop(terra::vect(x)) %>% 
      terra::mask(terra::vect(x)) %>%
      terra::xyFromCell(which(terra::values(.) == 0)) %>%
      as.data.frame() %>% 
      sf::st_as_sf(coords = c("x","y"), crs = sf::st_crs(x)) %>% 
      sf::st_set_geometry(sf_column)
    rm(x_bbox, sf_column)
  } else if (as.character(sf::st_geometry_type(x, by_geometry = FALSE)) == "MULTIPOINT") {
    observer <- sf::st_cast(x, "POINT")
  }
  
  return(x)
}