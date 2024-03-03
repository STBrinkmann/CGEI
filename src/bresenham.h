#ifndef BRESENHAM
#define BRESENHAM

#include <Rcpp.h>

extern int Sign2(const int dxy);
extern Rcpp::IntegerMatrix bresenham_map(const int x0, const int y0, const int radius, const int nc);
extern Rcpp::IntegerVector LoS_reference(const int x0_ref, const int y0_ref, const int r, const int nc_ref);
extern Rcpp::IntegerVector shared_LoS(const int r, const Rcpp::IntegerVector &los_ref_vec);

#endif