#' sf to raster interpolation
#' @description Function for computing a continuous raster map from a sf object by interpolating the sf values. For the interpolation Inverse Distance Weighting (IDW) is beeing used.
#'
#' @param observer object of class \code{sf}; Observer location(s) used for the interpolation. See ‘Details’ for valid sf geometry types
#' @param v character; Name of a variable in \code{observer} which values are used for the interpolation
#' @param aoi optional; object of class \code{sf}; Area of interest the raster map should me masked to. If NULL, the bounding box of \code{observer} will be used instead
#' @param max_distance numeric; Only observer locations within the buffer distance from the prediction location are used for interpolation 
#' @param n numeric; Number of nearest observer locations that should be used for the interpolation
#' @param beta numeric; Inverse istance power, β, determines the degree to which the nearer observer locations are preferred over more distant points
#' @param raster_res numeric; Resolution of the interpolatet raster map that is returned
#' @param spacing optional; numeric; Only if \code{observer} is a linestring (or polygon), points on the line (or on a grid) will be generated. The \code{spacing} parameter sets the distance in between the points on the line/grid
#' @param na_only logical; If TRUE, only NA values will be interpolated.
#' @param cores numeric; The number of cores to use, i.e. at most how many child processes will be run simultaneously
#' @param progress logical; Show progress bar?
#'
#' @details 
#' \code{observer} needs to be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON. If \code{observer} is a LINESTRING or MULTILINESTRING, 
#' points will be generated along the line(s) every \code{raster_res} meters. If \code{observer} is a POLYGON or MULTIPOLYGON, a grid with resolution = \code{raster_res} 
#' will be generated as the obersver locations.
#' The CRS (\code{\link[sf]{st_crs}}) needs to have a metric unit! 
#'
#' @return object of class \code{\link[terra]{rast}} 
#' @export
#' @importFrom sf st_crs st_geometry_type st_as_sf st_bbox
#' @importFrom dplyr rename
#' @importFrom raster raster
#' @importFrom terra rast vect values
sf_interpolat_IDW <- function(observer, v, aoi = NULL, max_distance = Inf,
                              n = Inf, beta = 2, raster_res = NULL,
                              spacing = raster_res, na_only = FALSE, 
                              cores = 1, progress = FALSE) {
  #### 1. Check input ####
  # Check observer
  valid_sf_types <- c("POINT", "MULTIPOINT", "LINESTRING", "MULTILINESTRING", "POLYGON", "MULTIPOLYGON")
  checkmate::assert(is(observer, "sf"), "observer must be an sf object")
  checkmate::assert(!sf::st_is_longlat(sf::st_crs(observer)), "observer must have a metric CRS")
  checkmate::assert(sf::st_geometry_type(observer, by_geometry = FALSE) %in% valid_sf_types, "observer must be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON")
  
  # v
  checkmate::assert(v %in% names(observer), "v must be a column name of observer")
  
  # aoi
  if(!is.null(aoi)){
    checkmate::assert(is(aoi, "sf"), "aoi must be a sf object")
    checkmate::assert(sf::st_geometry_type(aoi, by_geometry = FALSE) %in% c("POLYGON", "MULTIPOLYGON"), "observer must be POLYGON or MULTIPOLYGON")
    checkmate::assert(sf::st_crs(aoi) == sf::st_crs(observer), "aoi must have the same CRS as observer")
  }
  
  # max_distance
  checkmate::assert(max_distance > 0, "max_distance must be greater 0")
  mode = 1
  if (is.infinite(max_distance)) {
    mode = 0
  }
  
  # n
  checkmate::assert(n > 0, "n must be greater 0")
  
  # raster_res
  checkmate::assert(!is.null(raster_res) & raster_res > 0, "raster_res must be greater 0")
  
  # spacing
  if (is.null(spacing)) {
    spacing <- raster_res
  }
  
    #### 2. Convert observer to points ####
  if(progress) {
    message("Preprocessing:")
    pb = txtProgressBar(min = 0, max = 3, initial = 0, style = 3)
  }
  observer <- sf_to_POINT(observer, spacing, dsm_rast)
  
  if (progress) setTxtProgressBar(pb, 1)
  
  #### 3. Prepare data for interpolation analysis ####
  # observer
  if(is.null(aoi)){
    aoi <- sf::st_as_sfc(sf::st_bbox(observer))
  }
  
  # n
  if(is.infinite(n)){
    n = nrow(observer)
  }
  
  # mode
  if(mode == 1) {
    box_size <- (as.integer(max_distance/raster_res)) * (as.integer(max_distance/raster_res))
    mode <- ifelse(box_size > n, 1, 0)
  }
  
  # IWD raster
  iwd_rast <- terra::rast(ext = terra::ext(aoi),
                          crs = sf::st_crs(aoi)$proj4string,
                          resolution = raster_res)
  obs_cells <- terra::cellFromXY(iwd_rast, sf::st_coordinates(observer))
  
  # Remove observer outside aoi
  observer <- observer[which(!is.na(obs_cells)), ]
  obs_cells <- na.omit(obs_cells)
  
  iwd_rast[obs_cells] <- dplyr::pull(observer, v)
  
  if (progress) setTxtProgressBar(pb, 2)
  
  iwd_cpp_rast <- iwd_rast %>% terra::rast() %>% raster::raster()
  iwd_raster_vec <- terra::values(iwd_rast, mat = FALSE)
  
  if (progress) setTxtProgressBar(pb, 3)
  
  #### 4. IWD ####
  if (progress) {
    message("\nComputing IWD interpolation:")
  }
  
  iwd_vals <- IDW_cpp(rast = iwd_cpp_rast, x = iwd_raster_vec,
                      sf_x = sf::st_coordinates(observer)[,1],
                      sf_y = sf::st_coordinates(observer)[,2],
                      sf_z = dplyr::pull(observer, v), 
                      n = n, b = beta, radius = max_distance,
                      na_only = na_only, ncores = cores, display_progress = progress)
  
  iwd_rast[] <- iwd_vals
  
  iwd_rast <- iwd_rast %>% 
    terra::crop(terra::vect(aoi), mask = TRUE)
  
  return(iwd_rast)
}