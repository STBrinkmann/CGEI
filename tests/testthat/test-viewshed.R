# Viewshed engine (C++) used by vgvi(), vvi() and viewshed_list().

test_that("observer cells and returned cell numbers use R's 1-based numbering", {
  # non-square raster with a non-round origin, flat surface
  dsm <- mk_rast(matrix(0, 9, 13), res = 2, xmin = 512003, ymax = 5403101)
  rows <- c(1, 1, 9, 9, 5, 5, 1, 9)
  cols <- c(1, 13, 1, 13, 7, 1, 7, 13)
  vs <- cpp_vvi(dsm, rows, cols, h0 = rep(1.7, 8), radius = 6)
  xy <- terra::xyFromCell(dsm, terra::cellFromRowCol(dsm, rows, cols))
  for (k in seq_along(rows)) {
    own <- terra::cellFromXY(dsm, xy[k, , drop = FALSE])
    expect_true(own %in% vs[[k]]$visible_cells, info = paste("observer", k))
    # every visible cell lies within the radius (6 m = 3 cells) of the observer
    vxy <- terra::xyFromCell(dsm, vs[[k]]$visible_cells)
    expect_true(all(sqrt((vxy[, 1] - xy[k, 1])^2 + (vxy[, 2] - xy[k, 2])^2) <= 6 + 1e-9))
    # and the set equals the naive reference
    ref <- ref_viewshed(terra::as.matrix(dsm, wide = TRUE), rows[k], cols[k], 1.7, 3)
    expect_identical(vs[[k]]$visible_cells, ref_cells(ref$visible))
    expect_identical(vs[[k]]$viewshed, ref_cells(ref$seen))
  }
  # flat surface: the whole circle (29 cells for r = 3) is visible from the centre
  expect_length(vs[[5]]$visible_cells, 29)
})

test_that("hand-checked configurations", {
  flat <- matrix(0, 11, 11)
  dsm <- mk_rast(flat)
  # observer enclosed by a 10 m wall ring: only itself and the 8 wall cells
  wall <- flat
  wall[5:7, 5:7] <- 10
  wall[6, 6] <- 0
  vs <- cpp_vvi(mk_rast(wall), 6, 6, 1.7, 5)[[1]]
  expect_length(vs$visible_cells, 9)
  expect_setequal(vs$visible_cells, as.vector(outer((4:6) * 11, 5:7, `+`)))
  # a pole directly east of the observer hides the cells behind it on the axis
  pole <- flat
  pole[6, 7] <- 10
  vs <- cpp_vvi(mk_rast(pole), 6, 6, 1.7, 3)[[1]]
  cell <- function(r, c) (r - 1) * 11 + c
  expect_true(cell(6, 7) %in% vs$visible_cells)
  expect_false(cell(6, 8) %in% vs$visible_cells)
  expect_false(cell(6, 9) %in% vs$visible_cells)
  expect_true(cell(6, 5) %in% vs$visible_cells)   # west is free
  expect_true(cell(6, 3) %in% vs$visible_cells)
  # observer whose eye level is below the surface sees only its own cell
  vs <- cpp_vvi(dsm, 6, 6, -1, 5)[[1]]
  expect_identical(vs$visible_cells, as.integer(cell(6, 6)))
  expect_length(vs$viewshed, 81)  # the potential viewshed is still the full circle
})

test_that("orientation: rows run north to south, columns west to east", {
  m <- matrix(0, 21, 21)
  m[, 14] <- 30  # wall east of the observer (column 11)
  dsm <- mk_rast(m)
  vs <- cpp_vvi(dsm, 11, 11, 1.7, 9)[[1]]
  rc <- terra::rowColFromCell(dsm, vs$visible_cells)
  expect_true(all(rc[, 2] <= 14))         # nothing visible behind the wall
  expect_true(any(rc[, 2] == 2))          # far west is visible
  m2 <- t(m)                              # the same wall south of the observer
  vs2 <- cpp_vvi(mk_rast(m2), 11, 11, 1.7, 9)[[1]]
  rc2 <- terra::rowColFromCell(mk_rast(m2), vs2$visible_cells)
  expect_true(all(rc2[, 1] <= 14))
  expect_setequal(paste(rc2[, 1], rc2[, 2]), paste(rc[, 2], rc[, 1]))
})

test_that("viewshed equals the naive reference on random terrain", {
  scenes <- list(
    list(dsm = random_dsm(45, 60, seed = 1), r = 9),
    list(dsm = random_dsm(45, 60, na_frac = 0.05, seed = 2), r = 9),   # NA heights
    list(dsm = random_dsm(40, 13, n_blocks = 6, seed = 3), r = 9),     # narrower than 2r+1
    list(dsm = outer(1:30, 1:50, function(i, j) 5 * sin(i / 7) + 3 * cos(j / 5)), r = 11)
  )
  for (s in seq_along(scenes)) {
    sc <- scenes[[s]]
    set.seed(10 + s)
    obs <- data.frame(row = sample.int(nrow(sc$dsm), 25, TRUE), col = sample.int(ncol(sc$dsm), 25, TRUE))
    obs <- obs[!is.na(sc$dsm[cbind(obs$row, obs$col)]), ]
    h0 <- pmin(sc$dsm[cbind(obs$row, obs$col)], 3) + stats::runif(nrow(obs), 0, 3)
    ref <- lapply(seq_len(nrow(obs)), function(k) ref_viewshed(sc$dsm, obs$row[k], obs$col[k], h0[k], sc$r))
    for (early_stop in c(TRUE, FALSE)) {
      vs <- cpp_vvi(mk_rast(sc$dsm), obs$row, obs$col, h0, sc$r, early_stop = early_stop)
      expect_identical(lapply(vs, `[[`, "visible_cells"), lapply(ref, function(x) ref_cells(x$visible)),
                       info = paste("scene", s, "early_stop", early_stop))
      expect_identical(lapply(vs, `[[`, "viewshed"), lapply(ref, function(x) ref_cells(x$seen)),
                       info = paste("scene", s, "early_stop", early_stop))
    }
  }
})

test_that("NA heights neither block the view nor become visible", {
  m <- matrix(0, 15, 15)
  m[8, 9:10] <- NA                    # NA cells east of the observer
  m[8, 11] <- 5                       # a visible object behind them
  dsm <- mk_rast(m)
  vs <- cpp_vvi(dsm, 8, 8, 1.7, 6)[[1]]
  na_cells <- which(is.na(t(m)))
  expect_false(any(na_cells %in% vs$visible_cells))
  expect_false(any(na_cells %in% vs$viewshed))
  expect_true(((8 - 1) * 15 + 11) %in% vs$visible_cells)
})

test_that("observers outside the raster or with NA cell are skipped", {
  dsm <- mk_rast(matrix(0, 10, 10))
  vs <- cpp_vvi(dsm, c(NA, 0, 5, 11), c(5, 5, NA, 5), rep(1.7, 4), 3)
  for (k in 1:4) {
    expect_length(vs[[k]]$visible_cells, 0)
    expect_length(vs[[k]]$viewshed, 0)
  }
  v <- cpp_vgvi(dsm, dsm, c(NA, 0, 5), c(5, 5, NA), rep(1.7, 3), 3)
  expect_true(all(is.na(v)))
})
