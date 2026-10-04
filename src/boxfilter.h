// Separable box filters used by gavi() (focal means) and lacunarity() (box
// masses).
//
// Plain C++ (no Rcpp / R API) so that the functions can run in OpenMP worker
// threads and be compiled stand-alone for the sanitizer tests (dev/cpp-tests).
//
// All filters cost O(1) per cell, independent of the window size (the original
// implementations visited every cell of every window: O(w^2) per cell).
// Every output element is computed by exactly one thread in a fixed order, so
// results do not depend on the number of threads. NaN (R's NA) values are
// ignored and counted separately.

#ifndef CGEI_BOXFILTER_H
#define CGEI_BOXFILTER_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cgei {

inline int omp_thread_num() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

// Window sums along rows: for every row and every output column c in
// [0, out_ncol), H[row * out_ncol + c] = sum of the non-NA values in columns
// [c - left, c + right] clipped to [0, ncol); if Hc is not null it receives the
// number of non-NA values in that window. Sums of integer-valued data are exact.
inline void box_rows(const double *v, const int nrow, const int ncol, const int left,
                     const int right, const int out_ncol, double *H, int *Hc,
                     const int nthreads) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int row = 0; row < nrow; ++row) {
    const double *in = v + static_cast<std::size_t>(row) * ncol;
    double *out = H + static_cast<std::size_t>(row) * out_ncol;
    int *outc = Hc ? Hc + static_cast<std::size_t>(row) * out_ncol : nullptr;
    double s = 0.0;
    int n = 0;
    // window of column 0
    const int hi = std::min(right, ncol - 1);
    for (int k = 0; k <= hi; ++k) {
      const double x = in[k];
      if (!std::isnan(x)) {
        s += x;
        ++n;
      }
    }
    if (out_ncol > 0) {
      out[0] = s;
      if (outc) outc[0] = n;
    }
    for (int c = 1; c < out_ncol; ++c) {
      const int add = c + right;
      if (add < ncol) {
        const double x = in[add];
        if (!std::isnan(x)) {
          s += x;
          ++n;
        }
      }
      const int rem = c - left - 1;
      if (rem >= 0) {
        const double x = in[rem];
        if (!std::isnan(x)) {
          s -= x;
          --n;
        }
      }
      out[c] = s;
      if (outc) outc[c] = n;
    }
  }
}

// Window sums along columns of the row sums H (nrow x ncol): for every output
// row r in [0, out_nrow) the running sum over rows [r - up, r + down] clipped
// to [0, nrow) is handed to consume(row, c0, c1, sum, cnt) for the column
// segment [c0, c1) (sum[k] / cnt[k] belong to column c0 + k; cnt is null if Hc
// is null). Columns are processed in stripes, one stripe per task.
template <class Consume>
inline void box_cols(const double *H, const int *Hc, const int nrow, const int ncol,
                     const int up, const int down, const int out_nrow, const int nthreads,
                     Consume &&consume) {
  const int SW = 128;  // stripe width (columns)
  const int nstripes = (ncol + SW - 1) / SW;
  const int nt = std::max(1, nthreads);
  std::vector<double> acc_buf(static_cast<std::size_t>(nt) * SW);
  std::vector<int> cnt_buf(Hc ? static_cast<std::size_t>(nt) * SW : 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
#endif
  for (int st = 0; st < nstripes; ++st) {
    const int tid = omp_thread_num();
    const int c0 = st * SW;
    const int c1 = std::min(ncol, c0 + SW);
    const int w = c1 - c0;
    double *acc = &acc_buf[static_cast<std::size_t>(tid) * SW];
    int *cnt = Hc ? &cnt_buf[static_cast<std::size_t>(tid) * SW] : nullptr;
    std::fill(acc, acc + w, 0.0);
    if (cnt) std::fill(cnt, cnt + w, 0);
    auto add_row = [&](const int row, const int sign) {
      const double *h = H + static_cast<std::size_t>(row) * ncol + c0;
      if (sign > 0) {
        for (int k = 0; k < w; ++k) acc[k] += h[k];
      } else {
        for (int k = 0; k < w; ++k) acc[k] -= h[k];
      }
      if (cnt) {
        const int *hc = Hc + static_cast<std::size_t>(row) * ncol + c0;
        if (sign > 0) {
          for (int k = 0; k < w; ++k) cnt[k] += hc[k];
        } else {
          for (int k = 0; k < w; ++k) cnt[k] -= hc[k];
        }
      }
    };
    const int hi = std::min(down, nrow - 1);
    for (int k = 0; k <= hi; ++k) add_row(k, +1);
    if (out_nrow > 0) consume(0, c0, c1, acc, cnt);
    for (int row = 1; row < out_nrow; ++row) {
      const int add = row + down;
      if (add < nrow) add_row(add, +1);
      const int rem = row - up - 1;
      if (rem >= 0) add_row(rem, -1);
      consume(row, c0, c1, acc, cnt);
    }
  }
}

// Sliding maximum (or minimum) over windows [c, c + w - 1] of a sequence of
// length n (van Herk / Gil-Werman: 3 comparisons per element for any w).
// Writes n - w + 1 values to out. g and h are scratch buffers of length n.
template <bool Max>
inline void sliding_extreme(const double *in, const std::ptrdiff_t stride, const int n,
                            const int w, double *out, const std::ptrdiff_t out_stride,
                            double *g, double *h) {
  auto better = [](const double a, const double b) { return Max ? std::max(a, b) : std::min(a, b); };
  for (int i = 0; i < n; ++i) {
    const double x = in[i * stride];
    g[i] = (i % w == 0) ? x : better(g[i - 1], x);
  }
  for (int i = n - 1; i >= 0; --i) {
    const double x = in[i * stride];
    h[i] = (i % w == w - 1 || i == n - 1) ? x : better(h[i + 1], x);
  }
  for (int c = 0; c + w <= n; ++c) out[c * out_stride] = better(h[c], g[c + w - 1]);
}

// Box maximum (Max = true) or minimum of all w x w boxes anchored at their
// top-left cell (row, col), row in [0, nrow - w], col in [0, ncol - w].
// NaN values are ignored (a box without any value yields -inf / +inf).
// out has (nrow - w + 1) x (ncol - w + 1) elements, row-major.
template <bool Max>
inline void box_extreme(const double *v, const int nrow, const int ncol, const int w,
                        double *out, const int nthreads) {
  const int onr = nrow - w + 1, onc = ncol - w + 1;
  if (onr <= 0 || onc <= 0) return;
  const double fill = Max ? -std::numeric_limits<double>::infinity()
                          : std::numeric_limits<double>::infinity();
  const int nt = std::max(1, nthreads);
  std::vector<double> tmp(static_cast<std::size_t>(nrow) * onc);
  const int nmax = std::max(nrow, ncol);
  std::vector<double> scratch(static_cast<std::size_t>(nt) * 3 * nmax);

  // 1. rows: tmp[row, c] = extreme of v[row, c .. c + w - 1]
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
  for (int row = 0; row < nrow; ++row) {
    double *buf = &scratch[static_cast<std::size_t>(omp_thread_num()) * 3 * nmax];
    double *line = buf, *g = buf + nmax, *h = buf + 2 * nmax;
    const double *in = v + static_cast<std::size_t>(row) * ncol;
    for (int c = 0; c < ncol; ++c) line[c] = std::isnan(in[c]) ? fill : in[c];
    sliding_extreme<Max>(line, 1, ncol, w, tmp.data() + static_cast<std::size_t>(row) * onc, 1, g, h);
  }

  // 2. columns: out[row, c] = extreme of tmp[row .. row + w - 1, c]
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
  for (int c = 0; c < onc; ++c) {
    double *buf = &scratch[static_cast<std::size_t>(omp_thread_num()) * 3 * nmax];
    double *g = buf + nmax, *h = buf + 2 * nmax;
    sliding_extreme<Max>(tmp.data() + c, onc, nrow, w, out + c, onc, g, h);
  }
}

inline bool any_nan(const double *v, const std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isnan(v[i])) return true;
  }
  return false;
}

}  // namespace cgei

#endif
