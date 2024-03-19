#include <Rcpp.h>
#include "rsinfo.h"

#include <vector>
#include <set>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "eta_progress_bar.h"

using namespace Rcpp;

// Helper function to compute the focal weighted mean for each layer
std::vector<double> compute_focal_weighted_mean(const RasterInfo &x_ras, const NumericMatrix &x_mat,
                                                const NumericMatrix &lac, int row, int col, bool na_rm) {
  std::vector<double> weightedSums(x_mat.ncol(), 0.0); // Store weighted sums for each layer
  std::set<int> distinctWs; // To store distinct window sizes
  
  for (int l = 0; l < lac.nrow(); ++l) {
    int layer = lac(l, 0) - 1; // Adjusting layer index to 0-based
    if (layer < 0 || layer >= x_mat.ncol()) continue; // Skip invalid layers
    
    int w = lac(l, 1);
    int r = (w - 1) / 2;
    double weight = lac(l, 2);
    distinctWs.insert(w); // Add window size to set of distinct window sizes
    
    double sum = 0.0;
    int count = 0;
    
    // Expand the loop to include negative indices if na_rm = false, leading to NA if encountered
    for (int i = row - r; i <= row + r; ++i) {
      for (int j = col - r; j <= col + r; ++j) {
        // Skip cells outside the raster boundaries if na_rm is true
        if (na_rm && (i < 0 || i >= x_ras.nrow || j < 0 || j >= x_ras.ncol)) continue;
        
        // Assign NA and break if encountering out-of-bound indices with na_rm = false
        if (!na_rm && (i < 0 || i >= x_ras.nrow || j < 0 || j >= x_ras.ncol)) {
          weightedSums[layer] = NA_REAL;
          break;
        }
        
        double value = x_mat(j + i * x_ras.ncol, layer);
        if (NumericMatrix::is_na(value)) {
          if (na_rm) continue; // Skip NA values if na_rm is true
          else {
            weightedSums[layer] = NA_REAL; // Assign NA and break if encountering NA with na_rm = false
            break;
          }
        }
        
        sum += value;
        count++;
      }
      if (NumericMatrix::is_na(weightedSums[layer])) break; // If NA was assigned, no need to continue
    }
    
    if (!NumericMatrix::is_na(weightedSums[layer]) && count > 0) {
      weightedSums[layer] += (sum / count) * weight; // Add weighted focal mean to the sum for this layer
    }
  }
  
  // Normalize weighted sums by the number of distinct window sizes
  for (size_t i = 0; i < weightedSums.size(); ++i) {
    if (!NumericMatrix::is_na(weightedSums[i])) {
      weightedSums[i] /= distinctWs.size();
    }
  }
  
  return weightedSums;
}

// [[Rcpp::export]]
NumericMatrix focal_sum(S4 &x, const NumericMatrix &x_mat, const NumericMatrix &lac,
                        const bool na_rm = true, const int ncores = 1, const bool display_progress = false) {
  RasterInfo x_ras(x);
  NumericMatrix result(x_ras.nrow * x_ras.ncol, x_mat.ncol());
  
  // Progress bar
  ETAProgressBar pb_ETA;
  Progress pb(x_ras.ncell, display_progress, pb_ETA);
  
#ifdef _OPENMP
  omp_set_num_threads(ncores);
#pragma omp parallel for collapse(2)
#endif
  for (int row = 0; row < x_ras.nrow; ++row) {
    for (int col = 0; col < x_ras.ncol; ++col) {
      if ( !pb.is_aborted() ) {
        Progress::check_abort(); 
        pb.increment();
        
        std::vector<double> layer_means = compute_focal_weighted_mean(x_ras, x_mat, lac, row, col, na_rm);
        for (size_t layer = 0; layer < layer_means.size(); ++layer) {
          result(col + row * x_ras.ncol, layer) = layer_means[layer];
        }
      }
    }
  }
  
  return result;
}