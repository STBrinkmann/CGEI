#ifndef RSINFO
#define RSINFO

#include <Rcpp.h>

// Raster geometry passed from R as c(xmin, xmax, ymin, ymax, nrow, ncol)
// (see raster_geometry() in R/utils.R).
struct RasterInfo {
  double xmin, xmax, ymin, ymax, res;
  int nrow, ncol, ncell;

  explicit RasterInfo(const Rcpp::NumericVector &geom) {
    if (geom.size() != 6) {
      Rcpp::stop("Raster geometry must be c(xmin, xmax, ymin, ymax, nrow, ncol).");
    }
    xmin = geom[0];
    xmax = geom[1];
    ymin = geom[2];
    ymax = geom[3];
    nrow = static_cast<int>(geom[4]);
    ncol = static_cast<int>(geom[5]);
    if (nrow < 1 || ncol < 1 || !(xmax > xmin) || !(ymax > ymin)) {
      Rcpp::stop("Invalid raster geometry.");
    }
    if (static_cast<double>(nrow) * ncol > 2147483647.0) {
      Rcpp::stop("Rasters with more than 2^31 - 1 cells are not supported.");
    }
    ncell = nrow * ncol;
    res = (xmax - xmin) / ncol;
  }
};

#endif
