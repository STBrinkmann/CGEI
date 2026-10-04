# Compare two versions measured with bench_vgvi.R / bench_vvi.R (Markdown tables).
#
# Usage: Rscript benchmarks/compare.R <results_dir> <old_tag> <new_tag> [<old_label> <new_label>]
#   reads vgvi_<tag>.rds, vgvi_long_<tag>.rds (optional) and vvi_<tag>.rds
#   e.g.  Rscript benchmarks/compare.R benchmarks/results 040 041 0.4.0 0.4.1
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
tags <- args[2:3]
labels <- if (length(args) >= 5) args[4:5] else tags
rd <- function(f) if (file.exists(file.path(dir, f))) readRDS(file.path(dir, f)) else NULL
fmt <- function(x, d = 3) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d))
speedup <- function(old, new) ifelse(is.na(old) | is.na(new), "-", sprintf("**%.1fx**", old / new))
md_table <- function(df) {
  cat("| ", paste(names(df), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(df)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) cat("| ", paste(df[i, ], collapse = " | "), " |\n", sep = "")
  cat("\n")
}

for (kind in c("vgvi", "vgvi_long")) {
  a <- rd(sprintf("%s_%s.rds", kind, tags[1]))
  b <- rd(sprintf("%s_%s.rds", kind, tags[2]))
  if (is.null(a) || is.null(b)) next
  m <- merge(a$timing, b$timing, by = c("scene", "radius", "cores", "n_obs"), suffixes = c("_old", "_new"))
  m <- m[order(m$scene, m$radius, m$cores), ]
  cat(sprintf("### VGVI: C++ core (VGVI_cpp), milliseconds per observer (%d observers)\n\n", m$n_obs[1]))
  df <- data.frame(m$scene, m$radius, m$cores, fmt(m$ms_per_obs_old), fmt(m$ms_per_obs_new),
                   speedup(m$ms_per_obs_old, m$ms_per_obs_new))
  names(df) <- c("scene", "max_distance", "threads", labels, "speed-up")
  md_table(df)
  cat(sprintf("### VGVI: complete vgvi() call for %d observers, seconds\n\n", m$n_obs[1]))
  df <- data.frame(m$scene, m$radius, m$cores, fmt(m$e2e_s_old, 2), fmt(m$e2e_s_new, 2),
                   speedup(m$e2e_s_old, m$e2e_s_new))
  names(df) <- c("scene", "max_distance", "threads", labels, "speed-up")
  md_table(df)
  same <- vapply(names(a$values), function(k) identical(a$values[[k]], b$values[[k]]), logical(1))
  cat(sprintf("VGVI values identical in %d of %d configurations.\n\n", sum(same), length(same)))
}

a <- rd(sprintf("vvi_%s.rds", tags[1]))
b <- rd(sprintf("vvi_%s.rds", tags[2]))
if (!is.null(a) && !is.null(b)) {
  m <- merge(a$timing, b$timing, by = c("radius", "cores", "n_obs"), suffixes = c("_old", "_new"))
  m <- m[order(m$radius, m$cores), ]
  cat(sprintf("### VVI (city, %d observers), seconds\n\n", m$n_obs[1]))
  df <- data.frame(m$radius, m$cores,
                   fmt(m$vvi_s_old, 2), fmt(m$vvi_s_new, 2), speedup(m$vvi_s_old, m$vvi_s_new),
                   fmt(m$cumulative_s_old, 2), fmt(m$cumulative_s_new, 2),
                   speedup(m$cumulative_s_old, m$cumulative_s_new),
                   fmt(m$viewshed_s_old, 2), fmt(m$viewshed_s_new, 2),
                   speedup(m$viewshed_s_old, m$viewshed_s_new))
  names(df) <- c("max_distance", "threads",
                 paste("vvi()", labels), "speed-up",
                 paste("cumulative", labels), "speed-up ",
                 paste("viewshed", labels), " speed-up")
  md_table(df)
  same <- with(m, mean_vvi_old == mean_vvi_new & cvvi_old == cvvi_new & sum_n_views_old == sum_n_views_new)
  cat(sprintf("VVI results (mean VVI, cumulative VVI, sum of n_views) identical in %d of %d configurations.\n\n",
              sum(same), length(same)))
}
