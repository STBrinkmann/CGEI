#include <Rcpp.h>

#include <vector>

#include "los_geometry.h"

// Reference lines of sight (see los_geometry.h): 8 * r lines with r steps
// each, holding 0-based cell ids of a (2r+1) x (2r+1) reference grid centred
// at (x0_ref, y0_ref), NA padded. Exported for the tests.
// [[Rcpp::export]]
Rcpp::IntegerVector LoS_reference(const int x0_ref, const int y0_ref, const int r,
                                  const int nc_ref) {
  if (r < 0) Rcpp::stop("r must be >= 0.");
  const std::vector<int> los = cgei::los_reference(x0_ref, y0_ref, r, nc_ref);
  return Rcpp::IntegerVector(los.begin(), los.end());  // kNaInt == NA_INTEGER
}
