test_that("lacunarity works", {
  # Create a SpatRast as input
  mat_sample <- matrix(data = c(
     1,1,0,1,1,1,0,1,0,1,1,0,
     0,0,0,0,0,1,0,0,0,1,1,1,
     0,1,0,1,1,1,1,1,0,1,1,0,
     1,0,1,1,1,0,0,0,0,0,0,0,
     1,1,0,1,0,1,0,0,1,1,0,0,
     0,1,0,1,1,0,0,1,0,0,1,0,
     0,0,0,0,0,1,1,1,1,1,1,1,
     0,1,1,0,0,0,1,1,1,1,0,0,
     0,1,1,1,0,1,1,0,1,0,0,1,
     0,1,0,0,0,0,0,0,0,1,1,1,
     0,1,0,1,1,1,0,1,1,0,1,0,
     0,1,0,0,0,1,0,1,1,1,0,1
     ), nrow = 12, ncol = 12, byrow = TRUE
  )
  x <- terra::rast(mat_sample)

  lac_res <- CGEI::lacunarity(x)

  testthat::expect_equal(lac_res$r, c(3, 5, 7))
  testthat::expect_equal(lac_res[["ln(r)"]], log(lac_res$r))
  testthat::expect_equal(round(lac_res$Lac, 2), c(1.09, 1.02, 1.01))
  testthat::expect_equal(lac_res[["ln(Lac)"]], log(lac_res$Lac))
  # and equal to the naive reference
  testthat::expect_equal(lac_res$Lac, sapply(c(3, 5, 7), function(w) ref_lacunarity(mat_sample, w, TRUE)),
                         tolerance = 1e-12)
})

test_that("lacunarity: hand-computed examples", {
  # binary: four 3 x 3 boxes with masses 5, 4, 4, 5 -> E[S^2] / E[S]^2 = 20.5 / 20.25
  b <- matrix(c(1, 1, 0, 0,
                1, 1, 0, 0,
                0, 0, 1, 1,
                0, 0, 1, 1), 4, 4, byrow = TRUE)
  expect_equal(lacunarity(terra::rast(b), r_vec = 3)$Lac, 82 / 81)
  # continuous (3 x 4): two boxes with ranges 4 and 8 -> 1 + var / mean^2 = 1 + 8 / 36
  m <- matrix(c(1, 2, 3, 4,
                2, 3, 4, 5,
                3, 4, 5, 10), 3, 4, byrow = TRUE)
  expect_equal(lacunarity(terra::rast(m), r_vec = 3)$Lac, 11 / 9)
  # the same raster transposed (rows / columns swapped) gives the same boxes
  expect_equal(lacunarity(terra::rast(t(m)), r_vec = 3)$Lac, 11 / 9)
})

test_that("lacunarity equals the naive reference (NA, non-square, continuous)", {
  set.seed(9)
  mb <- matrix(stats::rbinom(27 * 35, 1, 0.45), 27, 35)
  mb[1:4, ] <- NA                 # NA rows at the border
  mb[10:12, 20:22] <- NA          # NA hole
  mc <- matrix(stats::runif(27 * 35), 27, 35)
  mc[sample.int(27 * 35, 50)] <- NA
  w <- c(3, 5, 9, 15)
  expect_equal(lacunarity(terra::rast(mb), r_vec = w)$Lac,
               sapply(w, function(k) ref_lacunarity(mb, k, TRUE)), tolerance = 1e-12)
  expect_equal(lacunarity(terra::rast(mc), r_vec = w)$Lac,
               sapply(w, function(k) ref_lacunarity(mc, k, FALSE)), tolerance = 1e-12)
  # a box whose new rim is NA keeps the mass of its valid core (bug in <= 0.3.1)
  m <- matrix(NA_real_, 10, 10)
  m[1:3, 1:3] <- 1
  m[2, 2] <- 0
  m[6:10, 6:10] <- rep_len(c(0, 1), 25)
  expect_equal(lacunarity(terra::rast(m), r_vec = c(3, 5, 7))$Lac,
               sapply(c(3, 5, 7), function(k) ref_lacunarity(m, k, TRUE)), tolerance = 1e-12)
})

test_that("lacunarity does not depend on the order of r_vec", {
  set.seed(4)
  x <- terra::rast(matrix(sample(0:1, 400, TRUE), 20, 20))
  a <- lacunarity(x, r_vec = c(3, 5, 9))
  b <- lacunarity(x, r_vec = c(9, 5, 3))
  expect_equal(b$r, c(9, 5, 3))
  expect_equal(a$Lac, rev(b$Lac))
})

test_that("box sizes larger than the raster are dropped with a warning", {
  x <- terra::rast(matrix(stats::rbinom(10 * 30, 1, 0.5), 10, 30))
  expect_warning(l <- lacunarity(x, r_vec = c(3, 9, 11, 21)), "larger than the raster")
  expect_equal(l$r, c(3, 9))
  expect_error(suppressWarnings(lacunarity(x, r_vec = 21)), "No valid box size")
})
