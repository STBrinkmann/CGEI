# Line-of-sight geometry shared by vgvi(), vvi() and viewshed_list().

test_that("line-of-sight geometry is identical to CGEI 0.3.1 (golden data)", {
  golden <- utils::read.csv(testdata_path("los_reference_golden.csv"))
  for (r in unique(golden$r)) {
    expect_identical(CGEI:::LoS_reference(r, r, r, 2L * r + 1L),
                     as.integer(golden$cell[golden$r == r]),
                     info = paste("r =", r))
  }
})

test_that("C++ geometry equals an independent R implementation", {
  for (r in 1:40) {
    expect_identical(CGEI:::LoS_reference(r, r, r, 2L * r + 1L),
                     ref_los_encode(ref_los_lines(r), r),
                     info = paste("r =", r))
  }
})

test_that("lines of sight are well formed and cover the whole circle", {
  for (r in c(2, 3, 5, 8, 13, 30)) {
    lines <- ref_los_lines(r)
    expect_length(lines, 8 * r)
    ok <- vapply(lines, function(L) {
      d2 <- L[, "dr"]^2 + L[, "dc"]^2
      nrow(L) > 0 &&
        all(d2 <= r^2) &&                               # inside the circle
        all(diff(d2) > 0) &&                            # moving away from the observer
        max(abs(L[1, ])) == 1 &&                        # starts next to the observer
        all(pmax(abs(diff(L[, "dr"])), abs(diff(L[, "dc"]))) == 1)  # 8-connected path
    }, logical(1))
    expect_true(all(ok), info = paste("r =", r))
    covered <- unique(do.call(rbind, lines))
    disc <- expand.grid(dr = -r:r, dc = -r:r)
    disc <- disc[disc$dr^2 + disc$dc^2 <= r^2 & (disc$dr != 0 | disc$dc != 0), ]
    expect_equal(nrow(covered), nrow(disc), info = paste("r =", r))
    expect_setequal(paste(covered[, "dr"], covered[, "dc"]), paste(disc$dr, disc$dc))
  }
})

test_that("lines of sight are symmetric under rotation by 90 degrees", {
  r <- 11
  lines <- ref_los_lines(r)
  key <- function(L) paste(L[, "dr"], L[, "dc"], collapse = ";")
  for (l in seq_len(2 * r)) {
    L <- lines[[l]]
    rotated <- cbind(dr = -L[, "dc"], dc = L[, "dr"])  # (dr, dc) -> (-dc, dr)
    expect_identical(key(lines[[l + 2 * r]]), key(rotated))
  }
})
