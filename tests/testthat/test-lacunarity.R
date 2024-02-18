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
})
