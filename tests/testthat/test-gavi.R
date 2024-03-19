test_that("GAVI works", {
  # Create a random raster with 10x10 cells
  r <- matrix(runif(10*10), nrow = 10, ncol = 10)
  r <- terra::rast(r)
  
  # Set up lacunarity with r = 3, lac = 1 for no weigths
  lac <- dplyr::tibble(
    name = names(r),
    i = 1,
    r = 3,
    `ln(r)` = log(r),
    Lac = 1,
    `ln(Lac)` = log(Lac)
  )
  gavi <- CGEI::gavi(r, lac)
  
  # Compare to focal mean with r = 3
  focal_mean <- terra::focal(r, w = 3, fun = "mean", na.rm = TRUE)
  focal_mean <- CGEI:::reclassify_jenks(focal_mean, 9)
  
  testthat::expect_equal(as.integer(gavi[]), as.integer(focal_mean[]), 
                         info = "GAVI should be equal to focal mean with regular rasters")
  
  # Create a random raster with 20x10 cells
  r <- matrix(runif(20*10), nrow = 20, ncol = 10)
  r <- terra::rast(r)
  
  # Set up lacunarity with r = 3, lac = 1 for no weigths
  lac <- dplyr::tibble(
    name = names(r),
    i = 1,
    r = 3,
    `ln(r)` = log(r),
    Lac = 1,
    `ln(Lac)` = log(Lac)
  )
  gavi <- CGEI::gavi(r, lac)
  
  # Compare to focal mean with r = 3
  focal_mean <- terra::focal(r, w = 3, fun = "mean", na.rm = TRUE)
  focal_mean <- CGEI:::reclassify_jenks(focal_mean, 9)
  
  testthat::expect_equal(as.integer(gavi[]), as.integer(focal_mean[]), 
                         info = "GAVI should be equal to focal mean with irregular rasters")
})