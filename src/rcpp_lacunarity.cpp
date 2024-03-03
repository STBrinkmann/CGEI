#include <Rcpp.h>
using namespace Rcpp;

// Helper function to compute the number of iterations of the main
// for-loop in the rcpp_lacunarity function.
int max_iterations (const IntegerVector w_vec,
                    const int mode,
                    const int mat_width,
                    const int mat_height) {
  int N_r, r, w;
  int iters = w_vec.size();
  for (int i = 0; i < w_vec.size(); i++) {
    w = w_vec[i];
    
    if (mode == 1) {
      r = w;
      N_r = (mat_width - r + 1) * (mat_height - r + 1);
    } else {
      r = (w-1)/2;
      N_r = (mat_width-(2*r))*(mat_height-(2*r));
    }
    iters += N_r;
  }
  
  return(iters);
}


// Calculate Lacunarity based on Plotnik ()
double lacunarity (NumericVector box_masses,
                   const int fun,
                   const int N_r) {
  double lac;
  if (box_masses.size() > 1) {
    if (fun == 1) {
      // 1. Max number of box values
      int max_value = max(box_masses);
      
      // 2. Frequency distribution n(S,r)
      IntegerVector n_S_r(max_value+1, 0);
      for (int j = 0; j < box_masses.size(); j++) {
        n_S_r[box_masses[j]] += 1;
      }
      
      // 3. Probability distribution Q(S,r)
      NumericVector Q_S_r(max_value+1, 0.0);
      for (int k = 0; k < Q_S_r.size(); k++) {
        Q_S_r[k] = n_S_r[k] / double(box_masses.size());
      }
      
      // 4. First and second moments of Q(S,r): S*Q(S,r) and S^2*Q(S,r)
      NumericVector first_moment(max_value+1, 0.0);
      NumericVector second_moment(max_value+1,0.0);
      for (int S = 0; S < Q_S_r.size(); S++) {
        first_moment[S] = S*Q_S_r[S];
        second_moment[S]= S*Q_S_r[S]*S;
      }
      double Z_1 = sum(first_moment);
      double Z_2 = sum(second_moment);
      
      // 5. Lacunarity
      lac = Z_2/(Z_1*Z_1);
    } else {
      lac = 1 + (( sd(box_masses)  *  sd(box_masses) ) /
        ( mean(box_masses)*mean(box_masses) ));
      
    }
  } else {
    lac = NA_REAL;
  }
  
  return(lac);
}


#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "eta_progress_bar.h"

// [[Rcpp::export]]
NumericVector rcpp_lacunarity(const NumericMatrix &mat,
                              const IntegerVector &r_vec,
                              const int fun,
                              const int ncores=1,
                              const bool display_progress=false) {
  // Basic raster information
  const int mat_width = mat.ncol();
  const int mat_height = mat.nrow();
  
  // Vector to store box masses of previous r
  NumericVector prev_box_mass_min(mat_width*mat_height, NA_REAL);
  NumericVector prev_box_mass_max(mat_width*mat_height, NA_REAL);
  NumericVector prev_box_mass_sum(mat_width*mat_height, NA_REAL);
  
  // Output
  NumericVector output(r_vec.size());
  
  // Progress bar
  ETAProgressBar pb_ETA;
  Progress pb(r_vec.size(), display_progress, pb_ETA);

  // Begin main loop
  for (int j = 0; j < r_vec.size(); j++) {
    const int r = r_vec[j];
    const int N_r = (mat_width - r + 1) * (mat_height - r + 1);
    
    // Gliding box algorithm
    // Init shared vector for parallel loop
    NumericVector box_masses(N_r);
    
#if defined(_OPENMP)
    omp_set_num_threads(ncores);
#pragma omp parallel for shared(box_masses, prev_box_mass_min, prev_box_mass_max, prev_box_mass_sum)
#endif
    for (int i = 0; i < N_r; i++) {
      if (!pb.is_aborted()) {
        Progress::check_abort();
        
        // Get x/y from i
        const int y = trunc(i / (mat_width - r + 1));
        const int x = i - (y * (mat_width - r + 1));
        const int cell = y * mat_width + x;
        
        // Pull previous box mass value
        const double prev_min = prev_box_mass_min[cell];
        const double prev_max = prev_box_mass_max[cell];
        const double prev_sum = prev_box_mass_sum[cell];
        
        double cell_min = NumericVector::is_na(prev_min) ? R_PosInf : prev_min;
        double cell_max = NumericVector::is_na(prev_max) ? R_NegInf : prev_max;
        double cell_sum = NumericVector::is_na(prev_sum) ? 0.0 : prev_sum;
        
        bool prev = true;
        if(NumericVector::is_na(prev_min) && NumericVector::is_na(prev_max) && NumericVector::is_na(prev_sum)) {
          prev = false;
        }
        
        // Get box mass (window of r*r):
        // If fun == 1, box-sum will be calculated, else box-range
        int n = 0;
        if(prev) {
          for (int bx = r_vec[j-1]; bx < r; bx++) {
            for (int by=0; by < r; by++) {
              const double cell_value = mat(x+bx,y+by);
              if (!NumericVector::is_na(cell_value)) {
                // Update min
                if (cell_value < cell_min)
                  cell_min = cell_value;
                
                // Update max
                if (cell_value > cell_max)
                  cell_max = cell_value;
                
                // Update sum
                cell_sum += cell_value;
                
                // Update n
                n += 1;
              }
            }
          }
          for (int bx = 0; bx < r_vec[j-1]; bx++) {
            for (int by=r_vec[j-1]; by < r; by++) {
              const double cell_value = mat(x+bx,y+by);
              if (!NumericVector::is_na(cell_value)) {
                // Update min
                if (cell_value < cell_min)
                  cell_min = cell_value;
                
                // Update max
                if (cell_value > cell_max)
                  cell_max = cell_value;
                
                // Update sum
                cell_sum += cell_value;
                
                // Update n
                n += 1;
              }
            }
          }
        } else {
          for (int bx = 0; bx < r; bx++) {
            for (int by=0; by < r; by++) {
              const double cell_value = mat(x+bx,y+by);
              if (!NumericVector::is_na(cell_value)) {
                // Update min
                if (cell_value < cell_min)
                  cell_min = cell_value;
                
                // Update max
                if (cell_value > cell_max)
                  cell_max = cell_value;
                
                // Update sum
                cell_sum += cell_value;
                
                // Update n
                n += 1;
              }
            }
          }
        }
        
        // Box statistic
        if (n > 0) {
          if (fun == 1) {
            box_masses[i] = cell_sum;
            prev_box_mass_sum[cell] = cell_sum;
          } else {
            box_masses[i] = (cell_max - cell_min);
            prev_box_mass_min[cell] = cell_min;
            prev_box_mass_max[cell] = cell_max;
          }
        } else {
          box_masses[i] = NA_REAL;
          prev_box_mass_sum[cell] = NA_REAL;
          prev_box_mass_min[cell] = NA_REAL;
          prev_box_mass_max[cell] = NA_REAL;
        }
      }
    }
    
    // Remove NA and compute Lacunarity based on
    NumericVector bm_narm = wrap(na_omit(box_masses));
    
    // Compute lacunarity
    double lac = lacunarity(bm_narm, fun, N_r);
    
    output[j] = lac;
    
    pb.increment();
  }
  
  return(output);
}