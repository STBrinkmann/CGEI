#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "natural_breaks.h"

using namespace Rcpp;

// Natural breaks of x into k classes, replacing
// classInt::classIntervals(style = "jenks" | "fisher") (see natural_breaks.h):
//  * "jenks":        identical to classInt (exact port of its R code)
//  * "fisher":       same optimal partition as classInt's Fortran routine,
//                    O(k n log n) instead of O(k n^2)
//  * "fisher_exact": exact port of the Fortran routine (tests only)
// x must only contain finite values; 2 <= k <= length(x). Returns k + 1 breaks.
// [[Rcpp::export]]
Rcpp::NumericVector jenks_breaks_cpp(const Rcpp::NumericVector &x, const int k,
                                     const std::string &style = "fisher") {
  std::vector<double> d(x.begin(), x.end());
  for (double v : d) {
    if (!std::isfinite(v)) Rcpp::stop("x must only contain finite values.");
  }
  std::sort(d.begin(), d.end());
  if (k < 2 || static_cast<std::size_t>(k) > d.size())
    Rcpp::stop("k must be between 2 and the number of values.");
  std::vector<double> brks;
  if (style == "jenks") {
    brks = cgei::jenks_classint(d, k);
  } else if (style == "fisher") {
    cgei::FisherDC dc(d, k);
    brks = cgei::fisher_breaks(d, dc.solve());
  } else if (style == "fisher_exact") {
    brks = cgei::fisher_breaks(d, cgei::fisher_exact(d, k));
  } else {
    Rcpp::stop("Unknown style.");
  }
  return Rcpp::NumericVector(brks.begin(), brks.end());
}
