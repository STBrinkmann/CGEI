#include <Rcpp.h>

#include <cmath>
#include <set>
#include <vector>

#include "boxfilter.h"
#include "parallel_progress.h"
#include "rsinfo.h"

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Focal step of gavi(): for every layer the sum over its `lac` rows of
//   focal mean in a w x w window (w = lac[, 2]) * Lac (lac[, 3]),
// divided by the number of distinct window sizes. lac[, 1] is the 1-based
// layer index (column of x_mat). x_mat holds one layer per column with cells in
// row-major order (terra::values(x, mat = TRUE)).
//
// na_rm = TRUE : mean of the non-NA cells of the window that lie inside the
//                raster (windows without any such cell are skipped).
// na_rm = FALSE: NA as soon as the window leaves the raster or contains NA.
//
// Box sums make every window size O(1) per cell (the original summed every
// window cell by cell: O(w^2) per cell and window). Sums, counts and the
// accumulation order are the same as in the original, so integer-valued
// rasters give bit-identical results.
// [[Rcpp::export]]
NumericMatrix focal_sum(const NumericVector &x, const NumericMatrix &x_mat, const NumericMatrix &lac,
                        const bool na_rm = true, const int ncores = 1,
                        const bool display_progress = false) {
  const RasterInfo x_ras(x);
  const int nrow = x_ras.nrow, ncol = x_ras.ncol;
  const std::size_t ncell = static_cast<std::size_t>(nrow) * ncol;
  const int nlayer = x_mat.ncol();
  if (static_cast<std::size_t>(x_mat.nrow()) != ncell)
    Rcpp::stop("x_mat must have one row per raster cell.");
  if (lac.ncol() < 3) Rcpp::stop("lac must have the columns i, r and Lac.");

  // Valid lac rows (layer inside x_mat) and number of distinct window sizes
  struct LacRow { int layer, radius; double weight; };
  std::vector<LacRow> rows;
  std::set<int> distinct_w;
  for (int l = 0; l < lac.nrow(); ++l) {
    const int layer = static_cast<int>(lac(l, 0) - 1);  // R (1-based) -> C++ (0-based)
    if (layer < 0 || layer >= nlayer) continue;
    const int w = static_cast<int>(lac(l, 1));
    distinct_w.insert(w);
    rows.push_back({layer, (w - 1) / 2, lac(l, 2)});
  }
  const double n_w = static_cast<double>(distinct_w.size());

  NumericMatrix result(static_cast<int>(ncell), nlayer);  // zero-initialised
  const int nthreads = cgei::resolve_threads(ncores);
  cgei::ParallelProgress progress(rows.size(), display_progress);

  std::vector<double> H(ncell);
  std::vector<int> Hc;
  std::vector<char> layer_has_na(nlayer, -1);

  for (const LacRow &lr : rows) {
    const double *v = &x_mat(0, lr.layer);
    if (layer_has_na[lr.layer] < 0) layer_has_na[lr.layer] = cgei::any_nan(v, ncell) ? 1 : 0;
    const bool has_na = layer_has_na[lr.layer] == 1;
    if (has_na) Hc.resize(ncell);
    const int rad = lr.radius;
    const int full = (2 * rad + 1) * (2 * rad + 1);
    const double weight = lr.weight;
    double *out = &result(0, lr.layer);
    const double na_real = NA_REAL;

    cgei::box_rows(v, nrow, ncol, rad, rad, ncol, H.data(), has_na ? Hc.data() : nullptr,
                   nthreads);
    cgei::box_cols(H.data(), has_na ? Hc.data() : nullptr, nrow, ncol, rad, rad, nrow, nthreads,
                   [&](const int row, const int c0, const int c1, const double *sum,
                       const int *cnt) {
      const int r_lo = std::max(0, row - rad), r_hi = std::min(nrow - 1, row + rad);
      const bool rows_inside = row - rad >= 0 && row + rad < nrow;
      double *o = out + static_cast<std::size_t>(row) * ncol;
      for (int col = c0; col < c1; ++col) {
        const int c_lo = std::max(0, col - rad), c_hi = std::min(ncol - 1, col + rad);
        const int count = cnt ? cnt[col - c0] : (r_hi - r_lo + 1) * (c_hi - c_lo + 1);
        if (na_rm) {
          if (count > 0) o[col] += (sum[col - c0] / count) * weight;
        } else {
          if (std::isnan(o[col])) continue;  // stays NA
          const bool inside = rows_inside && col - rad >= 0 && col + rad < ncol;
          if (!inside || count != full) {
            o[col] = na_real;
          } else {
            o[col] += (sum[col - c0] / count) * weight;
          }
        }
      }
    });
    progress.tick();
  }
  progress.finish();

  // Normalise by the number of distinct window sizes
  for (int layer = 0; layer < nlayer; ++layer) {
    double *out = &result(0, layer);
    for (std::size_t i = 0; i < ncell; ++i) {
      if (!std::isnan(out[i])) out[i] /= n_w;
    }
  }
  return result;
}
