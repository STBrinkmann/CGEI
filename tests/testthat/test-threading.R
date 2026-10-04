# Multi-threading: results must not depend on the number of threads or on the
# order / grouping of the observers (no shared or leaked per-thread state).

scene_inputs <- function() {
  s <- load_scene()
  valid <- !is.na(s$expected$cell)
  h0 <- terra::extract(s$dtm, cbind(s$observers$x, s$observers$y))[, 1] + 1.7
  list(s = s, rows = s$expected$row[valid], cols = s$expected$col[valid], h0 = h0[valid])
}

test_that("OpenMP information is available", {
  info <- CGEI:::cgei_openmp_info()
  expect_type(info$openmp, "logical")
  expect_gte(info$max_threads, 1)
})

test_that("VGVI / VVI results are identical for 1, 2 and 4 threads", {
  d <- scene_inputs()
  for (fun in 1:3) {
    a1 <- cpp_vgvi(d$s$dsm, d$s$greenspace, d$rows, d$cols, d$h0, 30, fun = fun, cores = 1)
    expect_identical(cpp_vgvi(d$s$dsm, d$s$greenspace, d$rows, d$cols, d$h0, 30, fun = fun, cores = 2), a1)
    expect_identical(cpp_vgvi(d$s$dsm, d$s$greenspace, d$rows, d$cols, d$h0, 30, fun = fun, cores = 4), a1)
  }
  v1 <- cpp_vvi(d$s$dsm, d$rows, d$cols, d$h0, 30, cores = 1)
  expect_identical(cpp_vvi(d$s$dsm, d$rows, d$cols, d$h0, 30, cores = 4), v1)
  geom <- CGEI:::raster_geometry(d$s$dsm)
  vals <- terra::values(d$s$dsm, mat = FALSE)
  c1 <- CGEI:::VVI_count_cpp(geom, vals, d$cols, d$rows, d$h0, 30, ncores = 1L)
  c4 <- CGEI:::VVI_count_cpp(geom, vals, d$cols, d$rows, d$h0, 30, ncores = 4L)
  expect_identical(c4, c1)
  expect_identical(c1$n_visible, lengths(lapply(v1, `[[`, "visible_cells")))
  expect_identical(c1$n_viewshed, lengths(lapply(v1, `[[`, "viewshed")))
  # per-cell counts (accumulated by all threads)
  e1 <- CGEI:::VVI_cells_cpp(geom, vals, d$cols, d$rows, d$h0, 30, ncores = 1L)
  expect_identical(CGEI:::VVI_cells_cpp(geom, vals, d$cols, d$rows, d$h0, 30, ncores = 2L), e1)
  expect_identical(CGEI:::VVI_cells_cpp(geom, vals, d$cols, d$rows, d$h0, 30, ncores = 4L), e1)
  many <- CGEI:::VVI_cells_cpp(geom, vals, rep(d$cols, 25), rep(d$rows, 25), rep(d$h0, 25), 30, ncores = 4L)
  expect_identical(many$visible_count, e1$visible_count * 25L)
  expect_identical(many$viewshed_count, e1$viewshed_count * 25L)
})

test_that("observers are independent: batch == single calls, any order, many per thread", {
  d <- scene_inputs()
  n <- length(d$rows)
  batch <- cpp_vgvi(d$s$dsm, d$s$greenspace, d$rows, d$cols, d$h0, 30, fun = 2, cores = 4)
  single <- vapply(seq_len(n), function(k) {
    cpp_vgvi(d$s$dsm, d$s$greenspace, d$rows[k], d$cols[k], d$h0[k], 30, fun = 2, cores = 1)
  }, numeric(1))
  expect_identical(batch, single)
  set.seed(1)
  perm <- sample.int(n)
  expect_identical(cpp_vgvi(d$s$dsm, d$s$greenspace, d$rows[perm], d$cols[perm], d$h0[perm], 30,
                            fun = 2, cores = 4), batch[perm])
  # every thread processes many observers in sequence
  many <- cpp_vgvi(d$s$dsm, d$s$greenspace, rep(d$rows, 20), rep(d$cols, 20), rep(d$h0, 20), 30,
                   fun = 2, cores = 4)
  expect_identical(many, rep(batch, 20))
  vv <- cpp_vvi(d$s$dsm, rep(d$rows, 10), rep(d$cols, 10), rep(d$h0, 10), 30, cores = 4)
  expect_identical(vv, rep(cpp_vvi(d$s$dsm, d$rows, d$cols, d$h0, 30, cores = 1), 10))
})

test_that("focal step and lacunarity are identical for 1, 2 and 4 threads", {
  set.seed(2)
  m1 <- matrix(stats::rbinom(130 * 170, 1, 0.4), 130, 170)
  m2 <- matrix(stats::runif(130 * 170), 130, 170)
  m2[sample.int(length(m2), 500)] <- NA
  x <- c(terra::rast(m1), terra::rast(m2))
  geom <- CGEI:::raster_geometry(x)
  xm <- terra::values(x, mat = TRUE) * 1.0
  lac <- cbind(i = c(1, 2, 1, 2), r = c(3, 9, 33, 65), Lac = c(1.1, 1.2, 1.3, 1.4))
  for (na_rm in c(TRUE, FALSE)) {
    f1 <- CGEI:::focal_sum(geom, xm, lac, na_rm, 1L)
    expect_identical(CGEI:::focal_sum(geom, xm, lac, na_rm, 2L), f1)
    expect_identical(CGEI:::focal_sum(geom, xm, lac, na_rm, 4L), f1)
  }
  w <- c(3L, 5L, 17L, 65L, 129L)
  for (layer in 1:2) {
    v <- xm[, layer]
    fun <- if (layer == 1) 1L else 0L
    l1 <- CGEI:::rcpp_lacunarity(geom, v, w, fun, 1L)
    expect_identical(CGEI:::rcpp_lacunarity(geom, v, w, fun, 2L), l1)
    expect_identical(CGEI:::rcpp_lacunarity(geom, v, w, fun, 4L), l1)
  }
})

test_that("R wrappers give identical results for different numbers of cores", {
  s <- load_scene()
  v1 <- suppressMessages(vgvi(s$observers, s$dsm, s$dtm, s$greenspace, max_distance = 30, cores = 1))
  v4 <- suppressMessages(vgvi(s$observers, s$dsm, s$dtm, s$greenspace, max_distance = 30, cores = 4))
  expect_identical(v4$VGVI, v1$VGVI)
  for (mode in c("VVI", "cumulative", "viewshed")) {
    w1 <- suppressMessages(vvi(s$observers, s$dsm, s$dtm, max_distance = 30, mode = mode, cores = 1))
    w4 <- suppressMessages(vvi(s$observers, s$dsm, s$dtm, max_distance = 30, mode = mode, cores = 4))
    if (mode == "viewshed") {
      expect_identical(terra::values(w4), terra::values(w1))
    } else if (mode == "VVI") {
      expect_identical(w4$VVI, w1$VVI)
    } else {
      expect_identical(w4, w1)
    }
  }
  x <- c(s$greenspace, s$dsm)
  l1 <- lacunarity(x, cores = 1)
  l4 <- lacunarity(x, cores = 4)
  expect_identical(l4, l1)
  set.seed(42)
  g1 <- gavi(x, l1, cores = 1)
  set.seed(42)
  g4 <- gavi(x, l1, cores = 4)
  expect_identical(terra::values(g4), terra::values(g1))
})

test_that("progress reporting from parallel loops works", {
  d <- scene_inputs()
  out <- utils::capture.output(
    v <- CGEI:::VGVI_cpp(CGEI:::raster_geometry(d$s$dsm), terra::values(d$s$dsm, mat = FALSE),
                         CGEI:::raster_geometry(d$s$greenspace),
                         as.numeric(terra::values(d$s$greenspace, mat = FALSE)),
                         d$cols, d$rows, d$h0, 30, 3L, 1, 6, ncores = 4L, display_progress = TRUE),
    type = "message")
  expect_identical(v, cpp_vgvi(d$s$dsm, d$s$greenspace, d$rows, d$cols, d$h0, 30, cores = 1))
})
