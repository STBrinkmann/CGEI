#include <Rcpp.h>
#include "rsinfo.h"
using namespace Rcpp;

inline std::vector<int> findBestIndices(std::vector<double> &d, const int N)
{   
  std::vector<int> indices(d.size());
  std::iota(indices.begin(), indices.end(), 0); // Fill with 0,1,2,...
  
  // Use partial_sort to sort the first N elements by their associated distances in 'd'
  std::partial_sort(indices.begin(), indices.begin() + N, indices.end(),
                    [&d](int i, int j) { return d[i] < d[j]; });
  
  indices.resize(N); // Resize vector to contain only the first N indices
  return indices;
}

inline double calc_idw(std::vector<double> &d, std::vector<double> &v, const double b){
  double numerator = 0.0;
  double denominator = 0.0;
  
  // Sum from i to n
  std::vector<double> weights(d.size());
  for(std::size_t i = 0; i < d.size(); ++i){
    weights[i] = 1 / std::pow(d[i], b);
    numerator += v[i] * weights[i];
    denominator += weights[i];
  }
  
  return numerator/denominator;
}


#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "eta_progress_bar.h"

// [[Rcpp::export]]
NumericVector IDW_cpp(S4 &rast, const NumericVector &x,
                      const NumericVector &sf_x, const NumericVector &sf_y, const NumericVector &sf_z,
                      const size_t n, const double b, const double radius,
                      const bool na_only=false, const int ncores=1, const bool display_progress=false)
{
  // Basic raster information
  RasterInfo rast_info(rast);
  
  // Output
  NumericVector out(x.size(), NA_REAL);
  
  // Progress bar
  ETAProgressBar pb_ETA;
  Progress pb(x.size(), display_progress, pb_ETA);
  
  // Main loop: Loop over all values of the raster x
#if defined(_OPENMP)
  omp_set_num_threads(ncores);
#pragma omp parallel for schedule(dynamic) shared(out)
#endif
  for(int j = 0; j < x.size(); ++j){
    if ( pb.increment() ) {
      if (na_only && !NumericVector::is_na(x[j])) {
        out[j] = x[j];
        continue;
      }
      // 1. Convert j to row/col and X/Y coordinates
      // row col from cell
      const int row_j = j / rast_info.ncol;
      const int col_j = j - (row_j * rast_info.ncol);
      
      // XY from cell
      const double y_j = rast_info.ymax - (row_j + 0.5) * rast_info.res;
      const double x_j = rast_info.xmin + (col_j + 0.5) * rast_info.res;
      
      
      // 2. Calculate distance to all cells and store their values
      // Distance (d) and value (z) vector
      std::vector<double> d, z;
      d.reserve(sf_x.size());
      z.reserve(sf_x.size());
      
      // Iterate over all cells that are within the radius
      for(int i = 0; i < sf_x.size(); ++i){
        // Distance
        const double dist = sqrt(pow(x_j - sf_x[i], 2) + pow(y_j - sf_y[i], 2));
        
        // If distance <= 0, use a small value (resolution / 4) instead
        if(dist <= radius) {
          d.push_back(std::max(dist, rast_info.res / 4));
          z.push_back(sf_z[i]);
        }
      }
      
      if(!d.empty()) {
        
        // 3. Sort by distance and select top n
        std::vector<int> idx = findBestIndices(d, std::min(n, d.size()));
        std::vector<double> z_top_n(idx.size()), d_top_n(idx.size());
        for(size_t t = 0; t < idx.size(); ++t) {
          z_top_n[t] = z[idx[t]];
          d_top_n[t] = d[idx[t]];
        }
        
        // 4. Compute IDW
        out[j] = calc_idw(d_top_n, z_top_n, b);
      }
        
    }
  }
  
  return out;
}