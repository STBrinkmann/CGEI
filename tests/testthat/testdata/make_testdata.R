# Generates the test data in this folder. Run from tests/testthat/testdata:
#   Rscript make_testdata.R
#
# Scene ("scene_*"): 50 x 70 cells of 2 m, non-square and with a non-round
# origin, so that row/column mix-ups and R (1-based) / C++ (0-based) offsets
# change results. Contents: sloped and undulating terrain, buildings of
# different heights, tree crowns (greenspace), grass (greenspace, DSM = DTM),
# roads, a hole of NA heights and a few NA greenspace cells.
# Grids are stored as ESRI ASCII grids (plain text, CRS is set by
# helper-data.R: EPSG:25832).
#
# greenspace_fine: a 1 m greenspace raster whose origin is shifted by 0.25 m
# and which does not cover the eastern part of the scene (cells outside count
# as not green). Its pattern differs from the 2 m greenspace.
#
# observers.csv: observer locations covering interior, edge and corner cells,
# a point exactly on a cell boundary, a point under a tree crown (eye level
# below the DSM) and a point inside the NA hole.
#
# Expected values (expected_scene.csv) are computed with the naive R reference
# implementation in ../helper-reference.R (not with the package's C++ code).
#
# los_reference_golden.csv: line-of-sight geometry produced by the original
# C++ implementation of CGEI 0.3.1 (LoS_reference()), used to make sure the
# rewritten geometry code is bit-identical.

suppressPackageStartupMessages(library(terra))
source(file.path("..", "helper-reference.R"))

nr <- 50; nc <- 70; res <- 2
xmin <- 512003; ymax <- 5403101
crs <- "EPSG:25832"

# Terrain (rounded to cm)
dtm <- outer(seq_len(nr), seq_len(nc),
             function(r, c) 100 + 0.04 * c - 0.03 * r + 0.8 * sin(c / 9) * cos(r / 11))
dtm <- round(dtm, 2)

# Buildings
obj <- matrix(0, nr, nc)
obj[8:15, 10:20] <- 12
obj[30:42, 15:22] <- 8
obj[5:9, 45:60] <- 20
obj[20:24, 40:44] <- 15.5
obj[35:47, 50:52] <- 6

# Tree crowns
canopy <- matrix(0, nr, nc)
add_tree <- function(m, r0, c0, rad, h) {
  for (i in max(1, r0 - rad):min(nr, r0 + rad)) {
    for (j in max(1, c0 - rad):min(nc, c0 + rad)) {
      d <- sqrt((i - r0)^2 + (j - c0)^2)
      if (d <= rad) m[i, j] <- max(m[i, j], round(h * sqrt(1 - (d / (rad + 0.5))^2), 2))
    }
  }
  m
}
for (t in list(c(25, 8, 3, 9), c(45, 5, 2, 7), c(14, 30, 3, 11), c(18, 33, 2, 8),
               c(28, 60, 3, 10), c(40, 64, 2, 6), c(3, 35, 2, 9), c(33, 35, 3, 12))) {
  canopy <- add_tree(canopy, t[1], t[2], t[3], t[4])
}
canopy[obj > 0] <- 0

dsm <- dtm + pmax(obj, canopy)
dsm[44:46, 30:33] <- NA  # hole without height information

grass <- matrix(FALSE, nr, nc)
grass[26:50, 1:12] <- TRUE
grass[12:20, 25:35] <- TRUE
green <- (canopy > 0 | grass) * 1
green[obj > 0] <- 0
green[3:4, 66:67] <- NA  # unknown greenspace (counts as not green)

to_rast <- function(m, res_, xmin_, ymax_) {
  r <- rast(nrows = nrow(m), ncols = ncol(m), xmin = xmin_, xmax = xmin_ + ncol(m) * res_,
            ymin = ymax_ - nrow(m) * res_, ymax = ymax_, crs = crs)
  values(r) <- as.vector(t(m))
  r
}
dsm_r <- to_rast(dsm, res, xmin, ymax)
dtm_r <- to_rast(dtm, res, xmin, ymax)
green_r <- to_rast(green, res, xmin, ymax)

# 1 m greenspace, shifted origin, covering only x < xmin + 120 m
fine_nc <- 120; fine_nr <- 104
fine_xmin <- xmin + 0.25; fine_ymax <- ymax + 1.75
fine <- matrix(0, fine_nr, fine_nc)
fy <- fine_ymax - (seq_len(fine_nr) - 0.5)       # cell centre y
fx <- fine_xmin + (seq_len(fine_nc) - 0.5)       # cell centre x
for (i in seq_len(fine_nr)) {
  for (j in seq_len(fine_nc)) {
    # green in a diagonal band and in a checkerboard park (differs from the 2 m data)
    band <- abs((fx[j] - xmin) - 1.3 * (ymax - fy[i]) - 10) < 6
    park <- (fx[j] - xmin) > 60 && (fx[j] - xmin) < 100 && (ymax - fy[i]) > 50 &&
      ((floor(fx[j]) + floor(fy[i])) %% 2 == 0)
    fine[i, j] <- as.numeric(band || park)
  }
}
fine_r <- to_rast(fine, 1, fine_xmin, fine_ymax)

writeRaster(dsm_r, "scene_dsm.asc", overwrite = TRUE, NAflag = -9999)
writeRaster(dtm_r, "scene_dtm.asc", overwrite = TRUE, NAflag = -9999)
writeRaster(green_r, "scene_greenspace.asc", overwrite = TRUE, NAflag = -9999)
writeRaster(fine_r, "scene_greenspace_fine.asc", overwrite = TRUE, NAflag = -9999)
unlink(list.files(".", pattern = "\\.(prj|aux\\.xml)$"))

# Observers
cx <- function(col) xmin + (col - 0.5) * res  # centre of column col (1-based)
cy <- function(row) ymax - (row - 0.5) * res  # centre of row row (1-based)
obs <- data.frame(
  label = c("street_1", "street_2", "street_3", "park", "next_to_wall", "between_buildings",
            "corner_top_left", "corner_top_right", "corner_bottom_left", "corner_bottom_right",
            "edge_west", "edge_east", "edge_north", "edge_south",
            "on_cell_boundary", "under_canopy", "in_na_hole"),
  x = c(cx(25), cx(38), cx(55), cx(26), cx(21), cx(42),
        cx(1), cx(70), cx(1), cx(70),
        cx(1), cx(70), cx(30), cx(48),
        xmin + 2 * 27, cx(30), cx(31)),
  y = c(cy(27), cy(12), cy(30), cy(19), cy(11), cy(27),
        cy(1), cy(1), cy(50), cy(50),
        cy(22), cy(36), cy(1), cy(50),
        cy(33), cy(14), cy(45))
)
write.csv(obs, "observers.csv", row.names = FALSE)

# Expected values from the R reference (max_distance = 30 m -> r = 15 cells)
radius <- 30; r_cells <- radius / res; eye <- 1.7
dsm_m <- as.matrix(dsm_r, wide = TRUE)
green_m <- as.matrix(green_r, wide = TRUE)
xy_dsm <- xyFromCell(dsm_r, seq_len(ncell(dsm_r)))
fine_on_dsm <- matrix(extract(fine_r, xy_dsm)[, 1], nr, nc, byrow = TRUE)
fine_on_dsm[is.na(fine_on_dsm)] <- 0

rows <- rowFromY(dsm_r, obs$y); cols <- colFromX(dsm_r, obs$x)
h0 <- extract(dtm_r, cbind(obs$x, obs$y))[, 1] + eye
expected <- do.call(rbind, lapply(seq_len(nrow(obs)), function(k) {
  valid <- !is.na(dsm_m[rows[k], cols[k]])
  if (!valid) {
    return(data.frame(label = obs$label[k], row = rows[k], col = cols[k], cell = NA, n_visible = NA,
                      n_viewshed = NA, vvi = NA, vgvi_none = NA, vgvi_exponential = NA,
                      vgvi_logit = NA, vgvi_fine_none = NA))
  }
  vs <- ref_viewshed(dsm_m, rows[k], cols[k], h0[k], r_cells)
  rg <- ref_rings(vs$visible, green_m, rows[k], cols[k], res)
  rg_fine <- ref_rings(vs$visible, fine_on_dsm, rows[k], cols[k], res)
  data.frame(
    label = obs$label[k], row = rows[k], col = cols[k],
    cell = (rows[k] - 1) * nc + cols[k],
    n_visible = sum(vs$visible), n_viewshed = sum(vs$seen),
    vvi = sum(vs$visible) / sum(vs$seen),
    vgvi_none = ref_vgvi_from_rings(rg, radius, "none"),
    vgvi_exponential = ref_vgvi_from_rings(rg, radius, "exponential", m = 1, b = 6),
    vgvi_logit = ref_vgvi_from_rings(rg, radius, "logit", m = 0.5, b = 8),
    vgvi_fine_none = ref_vgvi_from_rings(rg_fine, radius, "none")
  )
}))
write.csv(expected, "expected_scene.csv", row.names = FALSE)
print(expected, digits = 4)

# Golden line-of-sight geometry (requires CGEI <= 0.3.1 installed as "CGEI";
# kept as is otherwise).
if (requireNamespace("CGEI", quietly = TRUE) &&
    utils::packageVersion("CGEI") <= "0.3.1") {
  golden <- do.call(rbind, lapply(c(1, 2, 3, 4, 5, 7, 10, 16), function(r) {
    v <- CGEI:::LoS_reference(r, r, r, 2L * r + 1L)
    data.frame(r = r, line = rep(seq_len(8 * r) - 1, each = r), step = rep(seq_len(r) - 1, 8 * r),
               cell = v)
  }))
  write.csv(golden, "los_reference_golden.csv", row.names = FALSE, na = "NA")
}
