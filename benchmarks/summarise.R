# Summarise benchmark results as Markdown tables.
# Usage: Rscript benchmarks/summarise.R [results_dir]   (default benchmarks/results)
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[1] else file.path("benchmarks", "results")
rd <- function(f) if (file.exists(file.path(dir, f))) readRDS(file.path(dir, f)) else NULL
fmt <- function(x, d = 3) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d))
speedup <- function(old, new) ifelse(is.na(old) | is.na(new), "-", sprintf("**%.0fx**", old / new))
md_table <- function(df) {
  cat("| ", paste(names(df), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(df)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) cat("| ", paste(df[i, ], collapse = " | "), " |\n", sep = "")
  cat("\n")
}

vn <- rd("vgvi_new.rds")
vo <- rd("vgvi_old.rds")
if (!is.null(vn) && !is.null(vo)) {
  m <- merge(vo$timing, vn$timing, by = c("scene", "radius", "cores", "n_obs"), suffixes = c("_old", "_new"))
  m <- m[order(m$scene, m$radius, m$cores), ]
  cat("### VGVI: C++ core (VGVI_cpp), milliseconds per observer\n\n")
  md_table(data.frame(scene = m$scene, `max_distance` = m$radius, threads = m$cores,
                      `0.3.1` = fmt(m$ms_per_obs_old), `0.4.0` = fmt(m$ms_per_obs_new),
                      `speed-up` = speedup(m$ms_per_obs_old, m$ms_per_obs_new), check.names = FALSE))
  cat("### VGVI: complete vgvi() call for", m$n_obs[1], "observers, seconds\n\n")
  md_table(data.frame(scene = m$scene, `max_distance` = m$radius, threads = m$cores,
                      `0.3.1` = fmt(m$e2e_s_old, 2), `0.4.0` = fmt(m$e2e_s_new, 2),
                      `speed-up` = speedup(m$e2e_s_old, m$e2e_s_new), check.names = FALSE))
  cat("### VGVI: mean value of the index (results change because of the bug fixes)\n\n")
  md_table(unique(data.frame(scene = m$scene, `max_distance` = m$radius,
                             `0.3.1` = fmt(m$mean_vgvi_old, 4), `0.4.0` = fmt(m$mean_vgvi_new, 4),
                             check.names = FALSE)))
}

gn <- rd("gavi_new.rds")
go <- Filter(Negate(is.null), list(rd("gavi_old_500.rds"), rd("gavi_old_1000.rds")))
if (!is.null(gn)) {
  old_t <- if (length(go)) do.call(rbind, lapply(go, `[[`, "timing")) else NULL
  m <- if (!is.null(old_t)) {
    merge(old_t, gn$timing, by = c("size", "cores", "windows"), all.y = TRUE, suffixes = c("_old", "_new"))
  } else {
    gn$timing
  }
  m <- m[order(m$size, m$cores), ]
  cat("### GAVI: seconds (two-layer raster, default lacunarity box sizes)\n\n")
  get <- function(col) if (col %in% names(m)) m[[col]] else rep(NA_real_, nrow(m))
  md_table(data.frame(
    raster = paste0(m$size, " x ", m$size), threads = m$cores, `box sizes` = m$windows,
    `lacunarity 0.3.1` = fmt(get("lacunarity_s_old"), 2), `lacunarity 0.4.0` = fmt(get("lacunarity_s_new"), 2),
    `speed-up` = speedup(get("lacunarity_s_old"), get("lacunarity_s_new")),
    `focal step 0.3.1` = fmt(get("focal_s_old"), 2), `focal step 0.4.0` = fmt(get("focal_s_new"), 3),
    `speed-up ` = speedup(get("focal_s_old"), get("focal_s_new")),
    `gavi() 0.3.1` = fmt(get("gavi_s_old"), 2), `gavi() 0.4.0` = fmt(get("gavi_s_new"), 2),
    ` speed-up` = speedup(get("gavi_s_old"), get("gavi_s_new")), check.names = FALSE))

  # equality of the outputs (same inputs)
  for (o in go) {
    for (key in names(o$outputs)) {
      a <- o$outputs[[key]]
      b <- gn$outputs[[key]]
      if (is.null(b)) next
      same_lac <- identical(a$lac$Lac, b$lac$Lac)
      rel <- max(abs(a$focal - b$focal) / pmax(abs(a$focal), 1e-300), na.rm = TRUE)
      cat(sprintf("- %s (threads %s): lacunarity identical: %s; focal step: layer 1 identical: %s, max. relative difference %.1e\n",
                  sub(" .*", "", key), sub(".* ", "", key), same_lac,
                  identical(a$focal[, 1], b$focal[, 1]), rel))
    }
  }
  cat("\n")
}

wn <- rd("vvi_new.rds")
wo <- rd("vvi_old.rds")
if (!is.null(wn) && !is.null(wo)) {
  m <- merge(wo$timing, wn$timing, by = c("radius", "cores", "n_obs"), suffixes = c("_old", "_new"))
  m <- m[order(m$radius, m$cores), ]
  cat("### VVI (city, ", m$n_obs[1], " observers), seconds\n\n", sep = "")
  md_table(data.frame(`max_distance` = m$radius, threads = m$cores,
                      `vvi() 0.3.1` = fmt(m$vvi_s_old, 2), `vvi() 0.4.0` = fmt(m$vvi_s_new, 2),
                      `speed-up` = speedup(m$vvi_s_old, m$vvi_s_new),
                      `cumulative 0.3.1` = fmt(m$cumulative_s_old, 2), `cumulative 0.4.0` = fmt(m$cumulative_s_new, 2),
                      ` speed-up` = speedup(m$cumulative_s_old, m$cumulative_s_new), check.names = FALSE))
}
