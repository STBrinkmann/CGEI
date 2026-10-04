#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "boxfilter.h"
#include "parallel_progress.h"
#include "rsinfo.h"

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

// Lacunarity of a set of box masses (NA already removed).
// fun == 1 (binary data, Plotnick et al. 1993): Z2 / Z1^2 from the frequency
// distribution of the integer box sums.
// fun != 1 (continuous data, Hoechstetter et al. 2011): 1 + var / mean^2.
// Kept bit-identical to the original implementation: the same floating-point
// operations in the same order (sd() and mean() are evaluated once instead of
// twice each; minimum, maximum and the integer frequency table are exact and
// computed in parallel).
double lacunarity(const NumericVector &box_masses, const int fun, const int nthreads) {
  const R_xlen_t n = box_masses.size();
  if (n <= 1) return NA_REAL;
  if (fun != 1) {
    const double s = sd(box_masses);
    const double m = mean(box_masses);
    return 1 + ((s * s) / (m * m));
  }
  const double *x = box_masses.begin();
  double min_value = R_PosInf, max_value_d = R_NegInf;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static) \
    reduction(min : min_value) reduction(max : max_value_d)
#endif
  for (R_xlen_t j = 0; j < n; j++) {
    min_value = std::min(min_value, x[j]);
    max_value_d = std::max(max_value_d, x[j]);
  }
  if (min_value < 0 || max_value_d > 67108864.0) {
    // Not a 0/1 raster: the frequency table would be invalid (negative masses)
    // or huge; use the equivalent moment formula E[S^2] / E[S]^2 instead.
    long double z1 = 0, z2 = 0;
    for (R_xlen_t j = 0; j < n; ++j) {
      z1 += x[j];
      z2 += static_cast<long double>(x[j]) * x[j];
    }
    z1 /= n;
    z2 /= n;
    return static_cast<double>(z2 / (z1 * z1));
  }
  // 1. Max box mass
  const int max_value = static_cast<int>(max_value_d);
  // 2. Frequency distribution n(S, r) (per-thread tables if they are small)
  std::vector<int> n_S_r(max_value + 1, 0);
  const int nt = static_cast<double>(max_value + 1) * nthreads <= static_cast<double>(n) ? nthreads : 1;
  if (nt > 1) {
    std::vector<std::vector<int>> local(nt, std::vector<int>(max_value + 1, 0));
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
    {
      std::vector<int> &h = local[cgei::omp_thread_num()];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (R_xlen_t j = 0; j < n; j++) h[static_cast<R_xlen_t>(x[j])] += 1;
    }
    for (int t = 0; t < nt; ++t)
      for (int S = 0; S <= max_value; S++) n_S_r[S] += local[t][S];
  } else {
    for (R_xlen_t j = 0; j < n; j++) n_S_r[static_cast<R_xlen_t>(x[j])] += 1;
  }
  // 3. Probability distribution Q(S, r)
  // 4. First and second moments of Q(S, r): S * Q(S, r) and S^2 * Q(S, r)
  double Z_1 = 0.0, Z_2 = 0.0;
  for (int S = 0; S <= max_value; S++) {
    const double Q = n_S_r[S] / double(n);
    Z_1 += S * Q;
    Z_2 += S * Q * S;
  }
  // 5. Lacunarity
  return Z_2 / (Z_1 * Z_1);
}

}  // namespace

// Gliding-box lacunarity of a raster for every box size w in r_vec.
// Box mass: sum (fun == 1) or max - min (fun != 1) of the non-NA cells of the
// w x w box; boxes without any non-NA cell are ignored. Box sizes larger than
// the raster give NA.
//
// All box masses of one size are computed with separable sliding windows in
// O(1) per cell (the original grew every box incrementally: O(w^2 - w_prev^2)
// per box). This also fixes two bugs of the incremental scheme: an all-NA rim
// turned a valid box into NA and reset its accumulated mass, and the result
// depended on r_vec being sorted ascending.
// [[Rcpp::export]]
NumericVector rcpp_lacunarity(const Rcpp::NumericVector &x, const Rcpp::NumericVector &x_values,
                              const IntegerVector &r_vec, const int fun,
                              const int ncores = 1, const bool display_progress = false) {
  const RasterInfo x_ras(x);
  const int nrow = x_ras.nrow, ncol = x_ras.ncol;
  const std::size_t ncell = static_cast<std::size_t>(nrow) * ncol;
  if (static_cast<std::size_t>(x_values.size()) != ncell)
    Rcpp::stop("x_values does not match the dimensions of x.");

  const double *v = x_values.begin();
  const bool has_na = cgei::any_nan(v, ncell);
  const int nthreads = cgei::resolve_threads(ncores);
  const double na_real = NA_REAL;

  NumericVector output(r_vec.size());
  cgei::ParallelProgress progress(r_vec.size(), display_progress);
  std::vector<double> H, extreme;
  std::vector<int> Hc, count;

  for (R_xlen_t j = 0; j < r_vec.size(); j++) {
    const int w = r_vec[j];
    if (IntegerVector::is_na(w) || w < 1 || w > nrow || w > ncol) {
      output[j] = NA_REAL;
      progress.tick();
      continue;
    }
    const int onr = nrow - w + 1, onc = ncol - w + 1;
    const std::size_t N_r = static_cast<std::size_t>(onr) * onc;
    NumericVector box_masses = Rcpp::no_init(static_cast<R_xlen_t>(N_r));
    double *bm = box_masses.begin();

    // Box sums (fun == 1) and/or numbers of non-NA cells per box
    if (fun == 1 || has_na) {
      H.resize(static_cast<std::size_t>(nrow) * onc);
      if (has_na) Hc.resize(static_cast<std::size_t>(nrow) * onc);
      if (fun != 1) count.resize(N_r);
      cgei::box_rows(v, nrow, ncol, 0, w - 1, onc, H.data(), has_na ? Hc.data() : nullptr,
                     nthreads);
      cgei::box_cols(H.data(), has_na ? Hc.data() : nullptr, nrow, onc, 0, w - 1, onr,
                     nthreads,
                     [&](const int row, const int c0, const int c1, const double *sum,
                         const int *cnt) {
        const std::size_t base = static_cast<std::size_t>(row) * onc;
        for (int col = c0; col < c1; ++col) {
          const int n = cnt ? cnt[col - c0] : w * w;
          if (fun == 1) {
            bm[base + col] = n > 0 ? sum[col - c0] : na_real;
          } else {
            count[base + col] = n;
          }
        }
      });
    }

    // Box ranges (fun != 1)
    if (fun != 1) {
      extreme.resize(N_r);
      cgei::box_extreme<true>(v, nrow, ncol, w, bm, nthreads);
      cgei::box_extreme<false>(v, nrow, ncol, w, extreme.data(), nthreads);
      const double *ex = extreme.data();
      const int *cnt = has_na ? count.data() : nullptr;
      const long long n_box = static_cast<long long>(N_r);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
      for (long long i = 0; i < n_box; ++i) {
        bm[i] = (!cnt || cnt[i] > 0) ? bm[i] - ex[i] : na_real;
      }
    }

    // Remove NA (only boxes without any valid cell are NA) and compute lacunarity
    if (has_na) {
      NumericVector bm_narm = wrap(na_omit(box_masses));
      output[j] = lacunarity(bm_narm, fun, nthreads);
    } else {
      output[j] = lacunarity(box_masses, fun, nthreads);
    }
    progress.tick();
  }
  progress.finish();
  return output;
}

// Number of distinct non-NA values of x, counted only up to `limit` (stops
// early). Used by lacunarity() to detect binary rasters without computing all
// unique values of a continuous raster.
// [[Rcpp::export]]
int n_distinct_upto(const Rcpp::NumericVector &x, const int limit) {
  std::vector<double> seen;
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    const double v = x[i];
    if (std::isnan(v)) continue;
    if (std::find(seen.begin(), seen.end(), v) == seen.end()) {
      seen.push_back(v);
      if (static_cast<int>(seen.size()) >= limit) break;
    }
  }
  return static_cast<int>(seen.size());
}
