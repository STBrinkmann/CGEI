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

# 4 x 5 raster with the values 1..20 (row by row):
#    1  2  3  4  5
#    6  7  8  9 10
#   11 12 13 14 15
#   16 17 18 19 20
hand_mat <- matrix(1:20, 4, 5, byrow = TRUE)
call_focal <- function(layers, lac, na_rm = TRUE, cores = 1L) {
  x <- do.call(c, lapply(layers, terra::rast))
  CGEI:::focal_sum(CGEI:::raster_geometry(x), terra::values(x, mat = TRUE) * 1.0,
                   as.matrix(lac), na_rm, cores)
}
cell <- function(row, col) (row - 1) * 5 + col

test_that("focal means: hand-computed values (borders, NA, na_rm)", {
  lac <- data.frame(i = 1, r = 3, Lac = 1)
  f <- call_focal(list(hand_mat), lac)[, 1]
  expect_equal(f[cell(1, 1)], (1 + 2 + 6 + 7) / 4)                 # corner: 4 cells
  expect_equal(f[cell(1, 2)], (1 + 2 + 3 + 6 + 7 + 8) / 6)         # edge: 6 cells
  expect_equal(f[cell(2, 2)], 7)                                   # interior: centre value
  expect_equal(f[cell(3, 4)], 14)
  expect_equal(f[cell(4, 5)], (14 + 15 + 19 + 20) / 4)
  # na_rm = FALSE: windows leaving the raster are NA
  g <- call_focal(list(hand_mat), lac, na_rm = FALSE)[, 1]
  expect_true(all(is.na(g[c(cell(1, 1:5), cell(4, 1:5), cell(2:3, 1), cell(2:3, 5))])))
  expect_equal(g[c(cell(2, 2:4), cell(3, 2:4))], c(7, 8, 9, 12, 13, 14))
  # an NA cell is skipped (na_rm = TRUE) or propagates (na_rm = FALSE)
  m <- hand_mat
  m[2, 2] <- NA
  expect_equal(call_focal(list(m), lac)[cell(1, 1), 1], (1 + 2 + 6) / 3)
  expect_true(is.na(call_focal(list(m), lac, na_rm = FALSE)[cell(2, 3), 1]))
  expect_equal(call_focal(list(m), lac, na_rm = FALSE)[cell(3, 4), 1], 14)
})

test_that("focal step: weights, several window sizes and layer mapping", {
  layer2 <- hand_mat + 100
  # rows deliberately not ordered by layer
  lac <- data.frame(i = c(2, 1, 2, 1), r = c(3, 3, 5, 5), Lac = c(1, 2, 3, 0.5))
  f <- call_focal(list(hand_mat, layer2), lac)
  # cell (2, 3): 3 x 3 mean = 8, 5 x 5 window clipped to the whole raster = 10.5
  expect_equal(f[cell(2, 3), 1], (2 * 8 + 0.5 * 10.5) / 2)
  expect_equal(f[cell(2, 3), 2], (1 * 108 + 3 * 110.5) / 2)
  # layers do not influence each other
  f1 <- call_focal(list(hand_mat, layer2 * 0), lac)
  expect_identical(f1[, 1], f[, 1])
})

test_that("focal step equals the naive reference (NA, two layers, both na_rm)", {
  set.seed(3)
  l1 <- matrix(stats::rbinom(23 * 31, 1, 0.4), 23, 31)
  l2 <- matrix(stats::runif(23 * 31), 23, 31)
  l2[sample.int(23 * 31, 40)] <- NA
  lac <- data.frame(i = c(2, 1, 1, 2, 1), r = c(5, 3, 7, 3, 21), Lac = c(1.3, 1.1, 1.7, 2.0, 0.9))
  for (na_rm in c(TRUE, FALSE)) {
    f <- call_focal(list(l1, l2), lac, na_rm)
    ref <- ref_gavi_focal(list(l1, l2), lac, na_rm)
    expect_equal(f[, 1], as.vector(t(ref[[1]])), tolerance = 1e-12)
    expect_equal(f[, 2], as.vector(t(ref[[2]])), tolerance = 1e-12)
  }
})

test_that("natural breaks: hand example and classInt equivalence", {
  v <- c(1, 2, 3, 10, 11, 12, 20, 21, 22)
  expect_equal(CGEI:::jenks_breaks(v, 3, "fisher"), c(1, 6.5, 16, 22))
  expect_equal(CGEI:::jenks_breaks(v, 3, "jenks"), c(1, 3, 12, 22))
  # fewer distinct values than classes: every value is its own class
  expect_equal(CGEI:::jenks_breaks(c(0, 0, 1, 1, 1), 9), c(-0.5, 0.5, 1.5))
  expect_equal(CGEI:::jenks_breaks(rep(2, 5), 9), c(2, 2))
  skip_if_not_installed("classInt")
  ci <- function(d, style) {
    suppressWarnings(classInt::classIntervals(d, 9, style = style, warnLargeN = FALSE)$brks)
  }
  # within-class sum of squares of the classes defined by breaks
  sse <- function(d, brks) {
    cl <- findInterval(d, brks, rightmost.closed = TRUE, left.open = TRUE)
    cl[cl == 0] <- 1
    sum(tapply(d, cl, function(z) sum((z - mean(z))^2)))
  }
  set.seed(11)
  continuous <- list(stats::runif(300), c(stats::rnorm(200), stats::rnorm(200, 5)), stats::rexp(700))
  for (d in continuous) {
    expect_equal(CGEI:::jenks_breaks(d, 9, "fisher"), ci(d, "fisher"), tolerance = 1e-12)
    expect_equal(CGEI:::jenks_breaks(d, 9, "jenks"), ci(d, "jenks"), tolerance = 1e-12)
    expect_equal(CGEI:::jenks_breaks_cpp(d, 9L, "fisher"), CGEI:::jenks_breaks_cpp(d, 9L, "fisher_exact"))
  }
  # data with many repeated values: (near-)ties may be resolved differently,
  # but the classification is always (at least) as good as classInt's
  discrete <- list(round(stats::runif(500) * 30) / 30, sample(0:12, 400, TRUE) / 12)
  for (d in discrete) {
    for (style in c("fisher", "jenks")) {
      expect_lte(sse(d, CGEI:::jenks_breaks(d, 9, style)), sse(d, ci(d, style)) * (1 + 1e-12))
    }
  }
})

test_that("reclassify_jenks() works with few distinct values", {
  r <- terra::rast(matrix(rep(c(0, 0.5, 1), length.out = 100), 10, 10))
  rc <- CGEI:::reclassify_jenks(r, 9)
  expect_setequal(unique(terra::values(rc, mat = FALSE)), 1:3)
  expect_equal(terra::values(rc, mat = FALSE), match(terra::values(r, mat = FALSE), c(0, 0.5, 1)))
})

test_that("gavi() returns classes 1..9 and masks NA cells", {
  set.seed(5)
  m1 <- matrix(stats::rbinom(40 * 60, 1, 0.3), 40, 60)
  m2 <- matrix(stats::runif(40 * 60), 40, 60)
  m1[1:3, 1:3] <- NA
  x <- c(terra::rast(m1), terra::rast(m2))
  lac <- lacunarity(x, r_vec = c(3, 5, 9))
  g <- gavi(x, lac)
  v <- terra::values(g, mat = FALSE)
  expect_true(all(v[!is.na(v)] %in% 1:9))
  expect_true(all(is.na(v[c(1, 2, 61)])))
})
