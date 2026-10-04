test_that("VVI works", {
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

  # Calculate VVI
  vvi_results <- CGEI::vvi(observers, dsm_rast, dtm_rast)

  # 13 cells within 200 m (2 cells): 12 resp. 13 of them are visible
  testthat::expect_equal(round(vvi_results$VVI, 3), c(0.923, 1.000))
  xy <- sf::st_coordinates(observers)
  dsm_m <- terra::as.matrix(dsm_rast, wide = TRUE)
  for (k in 1:2) {
    ref <- ref_viewshed(dsm_m, terra::rowFromY(dsm_rast, xy[k, 2]), terra::colFromX(dsm_rast, xy[k, 1]), 1.7, 2)
    testthat::expect_equal(vvi_results$VVI[k], sum(ref$visible) / sum(ref$seen))
    testthat::expect_equal(vvi_results$n_visible_cells[k], sum(ref$visible))
  }
})

test_that("vvi() reproduces the reference on the test scene", {
  s <- load_scene()
  exp <- s$expected[!is.na(s$expected$cell), ]
  v <- suppressMessages(vvi(s$observers, s$dsm, s$dtm, max_distance = 30))
  expect_identical(v$label, exp$label)
  expect_equal(v$VVI, exp$vvi, tolerance = 1e-14)
  expect_identical(as.integer(v$n_visible_cells), as.integer(exp$n_visible))
  # by_row for point features uses the full cell lists: same result
  vb <- suppressMessages(vvi(s$observers, s$dsm, s$dtm, max_distance = 30, by_row = TRUE))
  vb <- vb[!is.na(vb$VVI), ]
  expect_equal(vb$VVI, exp$vvi, tolerance = 1e-14)
})

test_that("cumulative and viewshed modes agree with the reference", {
  s <- load_scene()
  exp <- s$expected[!is.na(s$expected$cell), ]
  dsm_m <- terra::as.matrix(s$dsm, wide = TRUE)
  h0 <- terra::extract(s$dtm, cbind(exp$row, exp$col) * 0 +
                         terra::xyFromCell(s$dsm, exp$cell))[, 1] + 1.7
  ref <- lapply(seq_len(nrow(exp)), function(k) ref_viewshed(dsm_m, exp$row[k], exp$col[k], h0[k], 15))
  vis <- lapply(ref, function(x) ref_cells(x$visible))
  seen <- lapply(ref, function(x) ref_cells(x$seen))

  cvvi <- suppressMessages(vvi(s$observers, s$dsm, s$dtm, max_distance = 30, mode = "cumulative"))
  expect_equal(cvvi, length(unique(unlist(vis))) / length(unique(unlist(seen))))

  vs <- suppressMessages(vvi(s$observers, s$dsm, s$dtm, max_distance = 30, mode = "viewshed"))
  # vvi() crops the DSM; compare by coordinates
  xy_vis <- terra::xyFromCell(s$dsm, unlist(vis))
  n_views <- terra::extract(vs$n_views, xy_vis)[, 1]
  counts <- table(unlist(vis))
  expect_equal(as.vector(n_views[match(as.integer(names(counts)), unlist(vis))]), as.vector(counts))
  expect_equal(sum(terra::values(vs$n_views), na.rm = TRUE), length(unlist(vis)))
})

test_that("viewshed_list() returns the visible cells at the right place", {
  s <- load_scene()
  exp <- s$expected[!is.na(s$expected$cell), ]
  obs <- s$observers[!is.na(s$expected$cell), ]
  vl <- suppressMessages(viewshed_list(obs, s$dsm, s$dtm, max_distance = 30))
  expect_length(vl, nrow(exp))
  for (k in seq_along(vl)) {
    xy <- cbind(obs$x[k], obs$y[k])
    expect_equal(terra::extract(vl[[k]], xy)[, 1], 1)       # observer cell is visible
    v <- terra::values(vl[[k]], mat = FALSE)
    expect_equal(sum(v == 1, na.rm = TRUE), exp$n_visible[k])
    expect_equal(sum(!is.na(v)), exp$n_viewshed[k])
  }
})
