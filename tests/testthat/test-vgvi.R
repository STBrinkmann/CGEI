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

  # Reference: naive R implementation on the full rasters
  xy <- sf::st_coordinates(observers)
  dsm_m <- terra::as.matrix(dsm_rast, wide = TRUE)
  green_m <- terra::as.matrix(greenspace_rast, wide = TRUE)
  ref <- sapply(1:2, function(k) {
    ref_vgvi(dsm_m, green_m, terra::rowFromY(dsm_rast, xy[k, 2]),
             terra::colFromX(dsm_rast, xy[k, 1]), 1.7, 200, 100)
  })
  testthat::expect_equal(vgvi_results$VGVI, ref, tolerance = 1e-12)
  testthat::expect_equal(round(vgvi_results$VGVI, 3), c(0.271, 0.5))
})

# Flat terrain, observer in the centre of an 11 x 11 raster of 1 m cells and
# max_distance = 3 m: all 29 cells of the circle are visible. Distance rings
# (distance rounded to whole metres, at least 1):
#   ring 1: observer + 4 neighbours (1 m) + 4 diagonals (1.41 m)  ->  9 cells
#   ring 2: 4 cells at 2 m + 8 cells at 2.24 m                    -> 12 cells
#   ring 3: 4 cells at 2.83 m + 4 cells at 3 m                    ->  8 cells
# With a single green cell, VGVI (mode "none" = mean over the 3 rings) is
#   1/24 for a green cell in ring 3, 1/36 in ring 2, 1/27 in ring 1, 0 outside.
probe_cases <- data.frame(
  dr = c(0, 0, -3, 3, -2, -2, 2, 2, 0, 2, -1, 1, 0, 0, 3, -2),
  dc = c(3, -3, 0, 0, -2, 2, -2, 2, 2, 0, 2, 1, 0, 4, 1, -3),
  expected = c(rep(1 / 24, 8), 1 / 36, 1 / 36, 1 / 36, 1 / 27, 1 / 27, 0, 0, 0)
)

test_that("single green cells are counted in the right ring (R/C++ cell offsets)", {
  dsm <- mk_rast(matrix(0, 11, 11), res = 1, xmin = 700.5, ymax = 3000.25)
  for (k in seq_len(nrow(probe_cases))) {
    p <- probe_cases[k, ]
    g <- matrix(0, 11, 11)
    g[6 + p$dr, 6 + p$dc] <- 1
    green <- mk_rast(g, res = 1, xmin = 700.5, ymax = 3000.25)
    rings <- cpp_rings(dsm, green, 6, 6, 1.7, 3)[[1]]
    expect_identical(rings$ring, 1:3)
    expect_identical(rings$n_visible, c(9L, 12L, 8L))
    ring_of_probe <- max(1, floor(sqrt(p$dr^2 + p$dc^2) + 0.5))
    expected_green <- if (p$expected > 0) as.numeric(1:3 == ring_of_probe) else c(0, 0, 0)
    expect_identical(rings$green, expected_green, info = paste("probe", p$dr, p$dc))
    expect_equal(cpp_vgvi(dsm, green, 6, 6, 1.7, 3), p$expected, tolerance = 1e-15,
                 info = paste("probe", p$dr, p$dc))
  }
})

test_that("vgvi() end-to-end: probes in all directions", {
  dsm <- mk_rast(matrix(0, 11, 11), res = 1, xmin = 700.5, ymax = 3000.25)
  obs <- mk_observers(dsm, 6, 6)
  for (k in seq_len(nrow(probe_cases))) {
    p <- probe_cases[k, ]
    g <- matrix(0, 11, 11)
    g[6 + p$dr, 6 + p$dc] <- 1
    res <- vgvi(obs, dsm, dsm, mk_rast(g, res = 1, xmin = 700.5, ymax = 3000.25),
                max_distance = 3, mode = "none")
    expect_equal(res$VGVI, p$expected, tolerance = 1e-15, info = paste("probe", p$dr, p$dc))
  }
})

test_that("fully green / fully grey views give 1 / 0 for every mode and resolution", {
  for (res in c(1, 2, 5)) {
    n <- 2 * round(60 / res) + 5
    dsm <- mk_rast(matrix(0, n, n), res = res)
    obs <- mk_observers(dsm, c((n + 1) / 2, 3, n - 2), c((n + 1) / 2, n - 2, 3))
    for (mode in c("none", "exponential", "logit")) {
      green <- vgvi(obs, dsm, dsm, mk_rast(matrix(1, n, n), res = res), max_distance = 60, mode = mode)
      grey <- vgvi(obs, dsm, dsm, mk_rast(matrix(0, n, n), res = res), max_distance = 60, mode = mode)
      expect_equal(green$VGVI, rep(1, 3), info = paste(res, mode))
      expect_equal(grey$VGVI, rep(0, 3), info = paste(res, mode))
    }
  }
})

test_that("decay weights are the integral of the decay function per 1 m ring", {
  # exponential decay with m = 1, b = 6: integral of 1 / (1 + 6x) = log(1 + 6x) / 6
  w <- ref_decay_weights(1:3, 3, "exponential", m = 1, b = 6)
  expect_equal(w, diff(log(1 + 6 * (0:3) / 3)) / 6, tolerance = 1e-5)
  # VGVI of a green cell in ring 3 is its share of the weights of the non-empty rings
  dsm <- mk_rast(matrix(0, 11, 11))
  g <- matrix(0, 11, 11)
  g[6, 9] <- 1
  v <- cpp_vgvi(dsm, mk_rast(g), 6, 6, 1.7, 3, fun = 2L, m = 1, b = 6)
  expect_equal(v, (1 / 8) * w[3] / sum(w), tolerance = 1e-12)
  wl <- ref_decay_weights(1:3, 3, "logit", m = 0.5, b = 8)
  vl <- cpp_vgvi(dsm, mk_rast(g), 6, 6, 1.7, 3, fun = 1L, m = 0.5, b = 8)
  expect_equal(vl, (1 / 8) * wl[3] / sum(wl), tolerance = 1e-12)
})

test_that("rings without visible cells are ignored", {
  # 5 m cells: rings 1..4, 6, 8, 9 ... contain no cell at all; a fully green
  # view must still give 1 (it gave ~0.9 before 0.4.0)
  dsm <- mk_rast(matrix(0, 41, 41), res = 5)
  g <- mk_rast(matrix(1, 41, 41), res = 5)
  rings <- cpp_rings(dsm, g, 21, 21, 1.7, 100)[[1]]
  expect_true(all(rings$n_visible > 0))
  expect_false(4 %in% rings$ring)
  expect_equal(cpp_vgvi(dsm, g, 21, 21, 1.7, 100, fun = 2L), 1)
})

test_that("vgvi() reproduces the reference on the test scene", {
  s <- load_scene()
  exp <- s$expected[!is.na(s$expected$cell), ]
  expect_message(
    v_none <- vgvi(s$observers, s$dsm, s$dtm, s$greenspace, max_distance = 30, mode = "none"),
    "1 point has been removed"
  )
  expect_identical(v_none$label, exp$label)
  expect_equal(v_none$VGVI, exp$vgvi_none, tolerance = 1e-12)
  v_exp <- suppressMessages(vgvi(s$observers, s$dsm, s$dtm, s$greenspace, max_distance = 30,
                                 mode = "exponential", m = 1, b = 6))
  expect_equal(v_exp$VGVI, exp$vgvi_exponential, tolerance = 1e-12)
  v_log <- suppressMessages(vgvi(s$observers, s$dsm, s$dtm, s$greenspace, max_distance = 30,
                                 mode = "logit", m = 0.5, b = 8))
  expect_equal(v_log$VGVI, exp$vgvi_logit, tolerance = 1e-12)
  # observer under a tree crown: only its own (green) cell
  expect_equal(v_none$VGVI[v_none$label == "under_canopy"], 1)
})

test_that("greenspace on a different grid (finer, shifted, smaller extent)", {
  s <- load_scene()
  exp <- s$expected[!is.na(s$expected$cell), ]
  v <- suppressMessages(vgvi(s$observers, s$dsm, s$dtm, s$greenspace_fine, max_distance = 30))
  expect_equal(v$VGVI, exp$vgvi_fine_none, tolerance = 1e-12)
})

test_that("early termination does not change results", {
  s <- load_scene()
  exp <- s$expected[!is.na(s$expected$cell), ]
  h0 <- terra::extract(s$dtm, cbind(s$observers$x, s$observers$y))[, 1] + 1.7
  valid <- !is.na(s$expected$cell)
  for (fun in 1:3) {
    a <- cpp_vgvi(s$dsm, s$greenspace, s$expected$row[valid], s$expected$col[valid], h0[valid], 30,
                  fun = fun, early_stop = TRUE)
    b <- cpp_vgvi(s$dsm, s$greenspace, s$expected$row[valid], s$expected$col[valid], h0[valid], 30,
                  fun = fun, early_stop = FALSE)
    expect_identical(a, b)
  }
})

test_that("invalid arguments are rejected", {
  dsm <- mk_rast(matrix(0, 11, 11))
  obs <- mk_observers(dsm, 6, 6)
  expect_error(vgvi(obs, dsm, dsm, dsm, max_distance = 3, cores = 0))
  expect_error(vgvi(obs, dsm, dsm, dsm, max_distance = 3, cores = 2.5))
  expect_error(vgvi(obs, dsm, dsm, dsm, max_distance = 0.3), "at least 1")
})
