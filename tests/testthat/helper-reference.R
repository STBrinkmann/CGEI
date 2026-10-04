# Independent, deliberately naive R reference implementations used by the tests.
#
# None of these functions call into the package's C++ code: the line-of-sight
# geometry, the visibility walk, the VGVI ring aggregation, the focal means and
# the lacunarity are all re-implemented here in plain R, following the
# algorithm definitions as literally as possible (no shared-prefix reuse, no
# early termination, no sliding windows). The C++ implementations are tested
# against these functions.
#
# Conventions (identical to terra): matrices are indexed [row, col] with row 1
# at the top (ymax) and col 1 at the left (xmin); R cell numbers are 1-based
# and row-major: cell = (row - 1) * ncol + col.

# -- Line-of-sight geometry ---------------------------------------------------

# Bresenham lines of sight for a radius of r cells, in the order used by the
# C++ code: 8 * r lines, each a matrix with columns dr (row offset) and dc
# (column offset) of the cells visited when walking away from the observer.
# Lines 0..r cover the first octant (towards (dr, dc) = (r, 0..r)); the other
# seven octants are obtained by rotation / reflection.
ref_los_lines <- function(r) {
  stopifnot(r >= 1)
  # First octant: walk Y (rows) one cell per step, X (cols) by Bresenham.
  octant <- lapply(0:r, function(i) {
    err <- r %/% 2
    X <- 0L
    Y <- 0L
    xs <- integer(0)
    ys <- integer(0)
    repeat {
      Y <- Y + 1L
      err <- err + i
      if (err >= r) {
        X <- X + as.integer(i > 0)
        err <- err - r
      }
      if (X^2 + Y^2 > r^2) break
      xs <- c(xs, X)
      ys <- c(ys, Y)
    }
    list(x = xs, y = ys)
  })

  lines <- vector("list", 8 * r)
  put <- function(idx, dr, dc) {
    lines[[idx + 1]] <<- cbind(dr = as.integer(dr), dc = as.integer(dc))
  }
  for (i in 0:r) {
    x <- octant[[i + 1]]$x
    y <- octant[[i + 1]]$y
    put(0 * r + i, y, x)
    put(2 * r + i, -x, y)
    put(4 * r + i, -y, -x)
    put(6 * r + i, x, -y)
    if (i != 0 && i != r) {
      put(2 * r - i, x, y)
      put(4 * r - i, -y, x)
      put(6 * r - i, -x, -y)
      put(8 * r - i, y, -x)
    }
  }
  lines
}

# Encode lines the same way as CGEI:::LoS_reference(): a vector of length
# 8 * r * r holding reference-grid cell ids (0-based, row-major in a
# (2r+1) x (2r+1) grid centred on the observer), NA padded.
ref_los_encode <- function(lines, r) {
  nc <- 2L * r + 1L
  out <- rep(NA_integer_, 8L * r * r)
  for (i in seq_along(lines)) {
    L <- lines[[i]]
    if (nrow(L) > 0) {
      out[(i - 1L) * r + seq_len(nrow(L))] <- as.integer((L[, "dr"] + r) * nc + (L[, "dc"] + r))
    }
  }
  out
}

# -- Viewshed -----------------------------------------------------------------

# Naive viewshed: every line is walked from the observer outwards with its own
# horizon (largest tangent so far). A cell is visible if its tangent is strictly
# larger than the horizon. Cells outside the raster or with NA height are
# skipped (they neither block nor are visible). Visibility requires the
# observer's eye level h0 to be above the DSM at the observer cell; the
# observer cell itself is always visible.
#
# Returns logical matrices `visible` and `seen` (cells reached by a line of
# sight that are inside the raster and not NA; used by VVI).
ref_viewshed <- function(dsm, row0, col0, h0, r, lines = ref_los_lines(r)) {
  nr <- nrow(dsm)
  nc <- ncol(dsm)
  visible <- matrix(FALSE, nr, nc)
  seen <- matrix(FALSE, nr, nc)
  visible[row0, col0] <- TRUE
  seen[row0, col0] <- TRUE
  above <- isTRUE(h0 > dsm[row0, col0])
  for (L in lines) {
    horizon <- -9999
    for (s in seq_len(nrow(L))) {
      rr <- row0 + L[s, 1]
      cc <- col0 + L[s, 2]
      if (rr < 1L || rr > nr || cc < 1L || cc > nc) next
      h <- dsm[rr, cc]
      if (is.na(h)) next
      seen[rr, cc] <- TRUE
      if (!above) next
      tangent <- (h - h0) / sqrt(L[s, 1]^2 + L[s, 2]^2)
      if (tangent > horizon) {
        horizon <- tangent
        visible[rr, cc] <- TRUE
      }
    }
  }
  list(visible = visible, seen = seen)
}

# R (1-based) cell numbers of TRUE entries of a [row, col] matrix.
ref_cells <- function(m) {
  idx <- which(m, arr.ind = TRUE)
  sort(as.integer((idx[, 1] - 1L) * ncol(m) + idx[, 2]))
}

# -- VGVI -----------------------------------------------------------------------

# Distance-ring histogram of the visible cells: ring = distance to the observer
# in metres, rounded half away from zero, at least 1. `green` is a matrix of
# greenspace values aligned with the DSM grid (NA counts as 0).
ref_rings <- function(visible, green, row0, col0, res) {
  idx <- which(visible, arr.ind = TRUE)
  d <- res * sqrt((idx[, 1] - row0)^2 + (idx[, 2] - col0)^2)
  ring <- pmax(1, floor(d + 0.5))
  g <- green[idx]
  g[is.na(g)] <- 0
  n_visible <- tapply(rep(1L, length(ring)), ring, sum)
  green_sum <- tapply(g, ring, sum)
  data.frame(
    ring = as.integer(names(n_visible)),
    n_visible = as.integer(n_visible),
    green = as.numeric(green_sum[names(n_visible)])
  )
}

# Decay weight of each ring: trapezoidal integral (300 intervals) of the decay
# function over the normalised distance interval [(ring - 1) / R, ring / R].
ref_decay_weights <- function(ring, radius, mode, m, b, n = 300) {
  f <- switch(mode,
    logit = function(x) 1 / (1 + exp(b * (x - m))),
    exponential = function(x) 1 / (1 + (b * x^m))
  )
  vapply(ring, function(k) {
    upper <- k / radius
    lower <- upper - 1 / radius
    h <- (upper - lower) / n
    y <- f(lower + (0:n) * h)
    h / 2 * (y[1] + y[n + 1]) + sum(h * y[2:n])
  }, numeric(1))
}

# VGVI from a ring histogram: (decay weighted) mean over the rings that contain
# at least one visible cell of the proportion of green visible cells.
ref_vgvi_from_rings <- function(rings, radius, mode = "none", m = 1, b = 6) {
  rings <- rings[rings$n_visible > 0, , drop = FALSE]
  raw <- rings$green / rings$n_visible
  if (mode == "none") return(mean(raw))
  w <- ref_decay_weights(rings$ring, radius, mode, m, b)
  sum(raw * w) / sum(w)
}

# Full VGVI for one observer on [row, col] matrices.
ref_vgvi <- function(dsm, green, row0, col0, h0, radius, res,
                     mode = "none", m = 1, b = 6) {
  r <- round(radius / res)
  vs <- ref_viewshed(dsm, row0, col0, h0, r)
  rings <- ref_rings(vs$visible, green, row0, col0, res)
  ref_vgvi_from_rings(rings, radius, mode, m, b)
}

# -- GAVI focal step -------------------------------------------------------------

# Focal mean over a (2r+1) x (2r+1) window. na_rm = TRUE: mean of the non-NA
# cells inside the raster (NA if there are none). na_rm = FALSE: NA if the
# window leaves the raster or contains NA.
ref_focal_mean <- function(mat, r, na_rm = TRUE) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  out <- matrix(NA_real_, nr, nc)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      rows <- (i - r):(i + r)
      cols <- (j - r):(j + r)
      inside <- all(rows >= 1 & rows <= nr) && all(cols >= 1 & cols <= nc)
      v <- mat[rows[rows >= 1 & rows <= nr], cols[cols >= 1 & cols <= nc]]
      if (na_rm) {
        v <- v[!is.na(v)]
        out[i, j] <- if (length(v)) sum(v) / length(v) else NA_real_
      } else {
        out[i, j] <- if (inside && !anyNA(v)) sum(v) / length(v) else NA_real_
      }
    }
  }
  out
}

# GAVI focal step: for every layer the sum over its lac rows of
# focal_mean(w) * Lac, divided by the number of distinct window sizes.
# `layers` is a list of [row, col] matrices, `lac` a data.frame with columns
# i (1-based layer), r (window size w) and Lac.
ref_gavi_focal <- function(layers, lac, na_rm = TRUE) {
  k <- length(unique(lac$r[lac$i >= 1 & lac$i <= length(layers)]))
  lapply(seq_along(layers), function(li) {
    acc <- matrix(0, nrow(layers[[li]]), ncol(layers[[li]]))
    for (l in which(lac$i == li)) {
      fm <- ref_focal_mean(layers[[li]], (lac$r[l] - 1) %/% 2, na_rm)
      if (na_rm) {
        acc <- acc + ifelse(is.na(fm), 0, fm * lac$Lac[l])
      } else {
        acc <- acc + fm * lac$Lac[l]
      }
    }
    acc / k
  })
}

# -- Lacunarity ----------------------------------------------------------------

# Gliding-box lacunarity for a box size w. Box mass = sum (binary data) or
# max - min (continuous data) of the non-NA cells; boxes without any non-NA cell
# are ignored.
ref_lacunarity <- function(mat, w, binary) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  if (w > nr || w > nc) return(NA_real_)
  masses <- numeric(0)
  for (i in seq_len(nr - w + 1)) {
    for (j in seq_len(nc - w + 1)) {
      v <- mat[i:(i + w - 1), j:(j + w - 1)]
      v <- v[!is.na(v)]
      if (!length(v)) next
      masses <- c(masses, if (binary) sum(v) else max(v) - min(v))
    }
  }
  if (length(masses) <= 1) return(NA_real_)
  if (binary) {
    mean(masses^2) / mean(masses)^2
  } else {
    1 + stats::var(masses) / mean(masses)^2
  }
}
