#' Lacunarity
#' @description Computes the Lacunarity of a Spatial Raster object or a list of Spatial Raster objects.
#'
#' @param x `SpatRaster` or list of `SpatRaster` or character; Raster image that can be loaded using \code{\link[terra]{rast}}.
#' x can also be character, referring to a folder. In that case all files in the folder with the ending ".tif" will be used.
#' @param r_vec integer vector (optional); Vector of box diameter values. Every `r` must be odd and >2. It is recommended
#' that no r value is greater than half of the shorter dimension of `x`. If `r_vec` is NULL (default), `r_vec` will be generated automatically.
#' @param r_max integer (optional); Maximum value of `r`, `r_vec` will be cut off for values >r.
#' @param plot logical; Should the summary Lacunarity figure, showing Lacunarity of all SpatRasters be printed?
#' @param plot_path folder path; If not NULL, a folder path to save Lacunarity plots (see details)
#' @param progress logical; Show progress bar?
#' @param ncores numeric; The number of cores to use.
#'
#' @details
#' \code{lacunarity} is based on the algorithm for binary images provided by \href{https://doi.org/10.1007/BF00125351}{Plotnick et al. 1993}.
#' \href{https://doi.org/10.1016/j.ecocom.2011.01.001}{Hoechstetter et al. 2011} further applyed this algorithm on continuous raster images.
#' If \code{x} is a binary raster (e.g. Land Use) Lacunarity is beening calculated using the Plotnicks algorithm.
#'
#' @return \code{\link[tibble]{tibble}} containing all Lacunarity values: \code{name} and \code{i} refer to the name and index of the raster input.
#' \code{x} is the box diameter (a 5X5 box has a diameter of 5 with a radius of 2). \code{lac} indicates the Lacunarity value. 
#' @export
#'
#' @examples
#' # Create a SpatRast as input
#' mat_sample <- matrix(data = c(
#'    1,1,0,1,1,1,0,1,0,1,1,0,
#'    0,0,0,0,0,1,0,0,0,1,1,1,
#'    0,1,0,1,1,1,1,1,0,1,1,0,
#'    1,0,1,1,1,0,0,0,0,0,0,0,
#'    1,1,0,1,0,1,0,0,1,1,0,0,
#'    0,1,0,1,1,0,0,1,0,0,1,0,
#'    0,0,0,0,0,1,1,1,1,1,1,1,
#'    0,1,1,0,0,0,1,1,1,1,0,0,
#'    0,1,1,1,0,1,1,0,1,0,0,1,
#'    0,1,0,0,0,0,0,0,0,1,1,1,
#'    0,1,0,1,1,1,0,1,1,0,1,0,
#'    0,1,0,0,0,1,0,1,1,1,0,1
#'    ), nrow = 12, ncol = 12, byrow = TRUE
#' )
#' x <- terra::rast(mat_sample)
#'
#' lacunarity(x)
#'
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom dplyr tibble bind_rows arrange
#' @importFrom terra unique as.matrix nlyr
#' @importFrom ggplot2 ggplot aes geom_line labs scale_x_continuous scale_y_continuous theme_minimal element_text theme ggsave
#' @importFrom checkmate testClass testCharacter assertClass assertIntegerish assertNumeric assertFlag assertNumeric
#' @importFrom magrittr %>%
#' @useDynLib CGEI, .registration = TRUE
lacunarity <- function(x, r_vec = NULL, r_max = NULL, plot = FALSE, plot_path = NULL, progress = FALSE, ncores = 1L) {
  
  # 1. Check input ----------------------------------------------------------
  # Check if x is a list and contains only SpatRaster objects
  if (checkmate::testClass(x, "list")) {
    if (!all(sapply(x, function(rr) checkmate::testClass(rr, "SpatRaster")))) {
      stop("all elements of x must be SpatRaster objects. Use the terra package for loading raster images.")
    }
    
    r_list <- list()
    for (i in seq_along(x)) {
      rr <- x[[i]]
      if (terra::nlyr(rr) > 1) {
        for (j in 1:terra::nlyr(rr)) {
          r_list[[length(r_list) + 1]] <- rr[[j]]
        }
      } else {
        r_list[[length(r_list) + 1]] <- rr
      }
    }
    
    x <- r_list
  } else if (checkmate::testCharacter(x)) {
    r_paths <- list.files(path = x, pattern = "\\.tif$", full.names = TRUE)
    
    assertLength(r_paths, lower = 1, .var.name = "r_paths", add = "No .tif files in folder path x.")
    
    x <- lapply(r_paths, terra::rast)
  } else if (checkmate::testClass(x, "SpatRaster")) {
    if (terra::nlyr(x) > 1) {
      r_list <- list()
      for (j in 1:terra::nlyr(x)) {
        r_list[[length(r_list) + 1]] <- x[[j]]
      }
      x <- r_list
    }
  } else {
    checkmate::assertClass(x, "SpatRaster", msg = "x must be a SpatRaster object. Use the terra package for loading raster images")
  }
  
  # Check r_vec
  if (!is.null(r_vec)) {
    r_vec[r_vec < 3] <- NA
    r_vec <- na.omit(r_vec)
    checkmate::assertIntegerish(r_vec, any.missing = FALSE, lower = 3)
    
    invalid_r <- !(r_vec %% 2)
    r_vec[invalid_r] <- 2 * round(r_vec[invalid_r] / 2) + 1
    
    if (length(r_vec) == 0) {
      stop("The provided r_vec has length 0. Maybe add values >2")
    }
  }
  
  # Check r_max
  if (!is.null(r_max)) {
    checkmate::assertNumeric(r_max, len = 1, any.missing = FALSE, lower = 1, add = "r_max must be numeric and greater than 0. Only the first value will be used.")
    r_max <- r_max[1]
  }
  
  # Check plot flag
  checkmate::assertFlag(plot, add = "plot must be logical (TRUE/FALSE)")
  
  # plot_path
  if(!is.null(plot_path)) {
    # Tempdir
    temp_path <- tempfile(paste0("Lacunarity_", gsub("-", "_", Sys.Date()), "__"), plot_path)
    dir.create(temp_path, recursive = TRUE)
  }
  
  # Check progress flag
  checkmate::assertFlag(progress, add = "progress must be logical (TRUE/FALSE)")
  
  # Check ncores
  checkmate::assertNumeric(ncores, lower = 1, any.missing = FALSE, add = "ncores must be numeric and at least 1")
  ncores <- ifelse(ncores <= 1, 1L, floor(ncores))
  
  
  # 2. Main loop ------------------------------------------------------------
  
  # Output tibble
  out <- dplyr::tibble(
    name = as.character(),
    i = as.integer(),
    r = as.integer(),
    "ln(r)" = as.integer(),
    Lac = as.numeric(),
    "ln(Lac)" = as.numeric()
  )
  
  # Convert x to list
  if(!is.list(x)) x <- list(x)
  
  for (i in seq_along(x)) {
    # Get current raster
    this_x <- x[[i]]
    if (progress) {
      cat(paste("Computing Lacunarity for", names(this_x)[1]))
      cat("\n")
    }
    
    # Calculate r vector from raster dimension
    if (is.null(r_vec)) {
      set_r_vec_null <- TRUE
      max_r <- 1
      while (2^(max_r+1)+1 < round(min(dim(this_x)[1:2])/2)) {
        max_r <- max_r + 1
      }
      
      r_vec <- c(2^(1:max_r)+1, round(min(dim(this_x)[1:2])/4)*2+1)
    }
    
    # r_max
    if (!is.null(r_max)) {
      r_vec <- r_vec[r_vec <= r_max]
    }
    
    r_vec <- unique(as.integer(r_vec))
    
    # Is r binary? (e.g. Greenspace raster)
    lac_fun <- as.integer(nrow(terra::unique(this_x)) <= 2)
    
    # Convert raster to matrix
    rast_mat <- terra::as.matrix(this_x, wide = TRUE)
    
    # Calculate Lacunarity for all w
    this_lac <- rcpp_lacunarity(mat = rast_mat,
                                r_vec = r_vec,
                                fun = lac_fun,
                                ncores = ncores,
                                display_progress = progress)
    
    this_out <- dplyr::tibble(
      name = rep(names(this_x)[1], length(this_lac)),
      i = i,
      r = r_vec,
      "ln(r)" = log(r_vec),
      Lac = this_lac,
      "ln(Lac)" = log(this_lac)
    )
    
    out <- dplyr::bind_rows(out, this_out)
    
    if(set_r_vec_null) r_vec <- NULL
    
    if(progress) cat("\n")
  }
  
  
  # Plot --------------------------------------------------------------------
  if (plot) {
    if (length(unique(out$name)) != length(unique(out$i))) {
      x_names <- as.character(out$i)
    } else {
      x_names <- as.character(out$name)
    }
    
    #max_x <- ceiling(max(out$`ln(r)`, na.rm = TRUE))
    #max_y <- round(max(out$`ln(Lac)`, na.rm = TRUE), 1)+0.1
    
    out_na_rm <- na.omit(out)
    max_x <- NA
    max_y <- NA
    
    for (i in 1:nrow(out_na_rm)) {
      r_sort <- sort(out_na_rm$`ln(r)`, decreasing = TRUE)
      lac_sort <- sort(out_na_rm$`ln(Lac)`, decreasing = TRUE)
      
      if(is.na(max_x) && is.finite(r_sort[i])){
        if(i > 1){
          warning("Plot may be false due to infinte values. Please double check with the output table")
        }
        max_x <- ceiling(r_sort[i])
      }
      
      if(is.na(max_y) && is.finite(lac_sort[i])){
        if(i > 1){
          warning("Plot may be false due to infinte values. Please double check with the output table")
        }
        max_y <- round(lac_sort[i], 1) + 0.1
      }
    }
    
    if(is.na(max_y) || is.na(max_x)) {
      warning("Plot can't be created. Please look into the output table")
      break
    }
    
    p <- ggplot2::ggplot(data = out,
                         mapping = ggplot2::aes(x = `ln(r)`,
                                                y = `ln(Lac)`,
                                                colour = x_names),
                         environment = environment()) +
      ggplot2::geom_line(lwd = 0.8) +
      ggplot2::labs(x = "ln(r)",
                    y = "ln(lacunarity)",
                    colour = "Legend") +
      ggplot2::scale_x_continuous(breaks = seq(0, max_x, 1),
                                  limits = c(1, max_x)) +
      ggplot2::scale_y_continuous(breaks = seq(0, max_y, 0.1),
                                  limits = c(0, max_y)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                     axis.text = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 12),
                     legend.title = ggplot2::element_text(size = 14))
    
    print(p)
  }
  
  
  # Save plot ---------------------------------------------------------------
  if (!is.null(plot_path)) {
    # Save out as CSV
    write.csv(out, file = file.path(temp_path, "Lacunarity.csv"), row.names = FALSE)
    
    # Save summary ggplot
    if (length(unique(out$name)) != length(unique(out$i))) {
      x_names <- as.character(out$i)
    } else {
      x_names <- as.character(out$name)
    }
    
    #max_x <- ceiling(max(out$`ln(r)`, na.rm = TRUE))
    #max_y <- round(max(out$`ln(Lac)`, na.rm = TRUE), 1)+0.1
    
    out_na_rm <- na.omit(out)
    max_x <- NA
    max_y <- NA
    
    for (i in 1:nrow(out_na_rm)) {
      r_sort <- sort(out_na_rm$`ln(r)`, decreasing = TRUE)
      lac_sort <- sort(out_na_rm$`ln(Lac)`, decreasing = TRUE)
      
      if(is.na(max_x) && is.finite(r_sort[i])){
        if(i > 1){
          warning("Plot may be false due to infinte values. Please double check with the output table")
        }
        max_x <- ceiling(r_sort[i])
      }
      
      if(is.na(max_y) && is.finite(lac_sort[i])){
        if(i > 1){
          warning("Plot may be false due to infinte values. Please double check with the output table")
        }
        max_y <- round(lac_sort[i], 1) + 0.1
      }
    }
    
    if(is.na(max_y) || is.na(max_x)) {
      warning("Plot can't be created. Please look into the output table")
      break
    }
    
    p <- ggplot2::ggplot(data = out,
                         mapping = ggplot2::aes(x = `ln(r)`,
                                                y = `ln(Lac)`,
                                                colour = x_names),
                         environment = environment()) +
      ggplot2::geom_line(lwd = 0.8) +
      ggplot2::labs(x = "ln(r)", y = "ln(lacunarity)",
                    colour = "Legend") +
      ggplot2::scale_x_continuous(breaks = seq(0, max_x, 1), limits = c(1, max_x)) +
      ggplot2::scale_y_continuous(breaks = seq(0, max_y, 0.1), limits = c(0, max_y)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                     axis.text = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 12),
                     legend.title = ggplot2::element_text(size = 14))
    
    ggplot2::ggsave(filename = file.path(temp_path, "Lacunarity_summary.svg"), 
                    plot = p, device = "svg", 
                    width = 7, height = 6.5, units = "in")
    ggplot2::ggsave(filename = file.path(temp_path, "Lacunarity_summary.png"), 
                    plot = p, device = "png", 
                    width = 7, height = 6.5, units = "in", dpi = 500)
    
    # Save unique plots
    if (length(x) > 1) {
      for (i in seq_along(x)) {
        this_out <- out[out$i == i, ]
        
        if (length(unique(this_out$name)) != length(unique(this_out$i))) {
          x_names <- as.character(this_out$i)
        } else {
          x_names <- as.character(this_out$name)
        }
        
        #max_x <- ceiling(max(this_out$`ln(r)`, na.rm = TRUE))
        #max_y <- round(max(this_out$`ln(Lac)`, na.rm = TRUE), 1)+0.1
        
        out_na_rm <- na.omit(this_out)
        max_x <- NA
        max_y <- NA
        
        for (i in 1:nrow(out_na_rm)) {
          r_sort <- sort(out_na_rm$`ln(r)`, decreasing = TRUE)
          lac_sort <- sort(out_na_rm$`ln(Lac)`, decreasing = TRUE)
          
          if(is.na(max_x) && is.finite(r_sort[i])){
            if(i > 1){
              warning("Plot may be false due to infinte values. Please double check with the output table")
            }
            max_x <- ceiling(r_sort[i])
          }
          
          if(is.na(max_y) && is.finite(lac_sort[i])){
            if(i > 1){
              warning("Plot may be false due to infinte values. Please double check with the output table")
            }
            max_y <- round(lac_sort[i], 1) + 0.1
          }
        }
        
        if(is.na(max_y) || is.na(max_x)) {
          warning("Plot can't be created. Please look into the output table")
          break
        }
        
        p <- ggplot2::ggplot(data = this_out,
                             mapping = ggplot2::aes(x = `ln(r)`,
                                                    y = `ln(Lac)`,
                                                    colour = x_names),
                             environment = environment()) +
          ggplot2::geom_line(lwd = 0.8) +
          ggplot2::labs(x = "ln(r)", y = "ln(lacunarity)",
                        colour = "Legend") +
          ggplot2::scale_x_continuous(breaks = seq(0, max_x, 1), limits = c(1, max_x)) +
          ggplot2::scale_y_continuous(breaks = seq(0, max_y, 0.1), limits = c(0, max_y)) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                         axis.text = ggplot2::element_text(size = 12),
                         legend.text = ggplot2::element_text(size = 12),
                         legend.title = ggplot2::element_text(size = 14))
        
        ggplot2::ggsave(filename = file.path(temp_path, paste0("Lacunarity_", i, ".svg")), 
                        plot = p, device = "svg", 
                        width = 7, height = 6.5, units = "in")
        ggplot2::ggsave(filename = file.path(temp_path, paste0("Lacunarity_", i, ".png")), 
                        plot = p, device = "png", 
                        width = 7, height = 6.5, units = "in", dpi = 500)
      } 
    }
  }
  return(dplyr::arrange(out, i))
}