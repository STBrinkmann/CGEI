test_that("VGVI works", {
  # Simulate observer locations as sf points
  observers <- sf::st_as_sf(data.frame(lon = c(5, 5.1), lat = c(52, 52.1)),
                            coords = c("lon", "lat"), crs = 4326)

  # Transform to UTM zone 31N for metric units
  observers <- sf::st_transform(observers, 32631)

  # Create synthetic rasters for DSM, DTM, and greenspace around observers
  # This is a simplified example assuming flat terrain and random greenspaces
  bbox_observers <- sf::st_bbox(observers)
  x_range <- bbox_observers[c("xmin", "xmax")]
  y_range <- bbox_observers[c("ymin", "ymax")]

  # Create a raster with 1km x 1km around the observers with a resolution of 100m for DSM, DTM, and greenspace
  set.seed(123)
  dsm_rast <- terra::rast(res = 100,
                          xmin=min(x_range)-1000, xmax=max(x_range) + 1000,
                          ymin=min(y_range)-1000, ymax=max(y_range) + 1000,
                          crs = terra::crs(observers))
  dsm_rast[] <- runif(terra::ncell(dsm_rast), 0, 1.5) # Assign random heights

  dtm_rast <- terra::rast(dsm_rast, vals=0) # Flat terrain

  greenspace_rast <- terra::rast(dsm_rast)
  greenspace_rast[] <- sample(0:1, terra::ncell(dsm_rast), replace=TRUE)

  # Calculate VGVI
  vgvi_results <- CGEI::vgvi(observers, dsm_rast, dtm_rast, greenspace_rast)
  
  testthat::expect_equal(round(vgvi_results$VGVI, 3), c(0.004, 0.001))
})
