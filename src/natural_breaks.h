// Natural breaks (Jenks / Fisher) used by reclassify_jenks() via jenks_breaks_cpp().
//
// Plain C++ (no Rcpp / R API) so that it can be compiled stand-alone for the
// sanitizer tests in dev/cpp-tests.
//
//  * jenks_classint(): exact port of classInt's style = "jenks" (pure R):
//    same arithmetic, same tie rule, same breaks (data values at the upper
//    class limits). O(n^2 k), but in C++.
//  * fisher_exact(): exact port of classInt's Fortran routine FISH (style =
//    "fisher"), O(k n^2). Only used to validate FisherDC.
//  * FisherDC: the same optimal partition (minimal within-class sum of squared
//    deviations, ties resolved towards the smallest class start) computed with
//    the divide-and-conquer optimisation of the dynamic programme:
//    O(k n log n) instead of O(k n^2).
//  * fisher_breaks(): classInt's breaks for a partition: minimum, midpoints
//    between adjacent classes, maximum.
//
// Input vectors must be sorted ascending and contain only finite values.

#ifndef CGEI_NATURAL_BREAKS_H
#define CGEI_NATURAL_BREAKS_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cgei {

// classInt style = "jenks" (R code ported 1:1, including its initialisation
// of the first row with 0 and the remaining rows with .Machine$double.xmax).
inline std::vector<double> jenks_classint(const std::vector<double> &d, const int k) {
  const int n = static_cast<int>(d.size());
  std::vector<int> mat1(static_cast<std::size_t>(n) * k, 1);
  std::vector<double> mat2(static_cast<std::size_t>(n) * k, 0.0);
  for (std::size_t i = k; i < mat2.size(); ++i) mat2[i] = DBL_MAX;
  auto at = [k](const int l, const int j) {  // 1-based (l, j) -> index
    return static_cast<std::size_t>(l - 1) * k + (j - 1);
  };
  double v = 0;
  for (int l = 2; l <= n; ++l) {
    double s1 = 0, s2 = 0, w = 0;
    for (int m = 1; m <= l; ++m) {
      const int i3 = l - m + 1;
      const double val = d[i3 - 1];
      s2 = s2 + val * val;
      s1 = s1 + val;
      w = w + 1;
      v = s2 - (s1 * s1) / w;
      const int i4 = i3 - 1;
      if (i4 != 0) {
        for (int j = 2; j <= k; ++j) {
          const double cand = v + mat2[at(i4, j - 1)];
          if (mat2[at(l, j)] >= cand) {
            mat1[at(l, j)] = i3;
            mat2[at(l, j)] = cand;
          }
        }
      }
    }
    mat1[at(l, 1)] = 1;
    mat2[at(l, 1)] = v;
  }
  // kclass: upper index of every class
  std::vector<int> kclass(k);
  for (int j = 0; j < k; ++j) kclass[j] = j + 1;
  kclass[k - 1] = n;
  int row = n;
  for (int j = k; j >= 1; --j) {
    if (row < 1) throw std::runtime_error("Degenerate natural breaks partition.");
    const int id = mat1[at(row, j)] - 1;
    if (j >= 2) kclass[j - 2] = id;
    row = id;
  }
  std::vector<double> brks;
  brks.push_back(d[0]);
  for (int j = 0; j < k; ++j) {
    if (kclass[j] < 1) throw std::runtime_error("Degenerate natural breaks partition.");
    brks.push_back(d[kclass[j] - 1]);
  }
  return brks;
}

// Breaks from the first index (1-based) of every class of a sorted vector:
// overall minimum, midpoints between adjacent classes, overall maximum.
inline std::vector<double> fisher_breaks(const std::vector<double> &d, const std::vector<int> &first) {
  const int k = static_cast<int>(first.size());
  const int n = static_cast<int>(d.size());
  std::vector<double> brks;
  brks.push_back(d[first[0] - 1]);
  for (int c = 0; c + 1 < k; ++c) {
    const double max_c = d[first[c + 1] - 2];   // last element of class c
    const double min_next = d[first[c + 1] - 1];
    brks.push_back((max_c + min_next) / 2);
  }
  brks.push_back(d[n - 1]);
  return brks;
}

// classInt's Fortran routine FISH, ported 1:1 (O(k n^2)).
inline std::vector<int> fisher_exact(const std::vector<double> &x, const int k) {
  const int m = static_cast<int>(x.size());
  const double big = static_cast<double>(10.E30f);  // R1MACH2 (single precision literal)
  std::vector<double> work(static_cast<std::size_t>(m) * k, big);
  std::vector<int> iwork(static_cast<std::size_t>(m) * k, 0);
  auto at = [k](const int i, const int j) { return static_cast<std::size_t>(i - 1) * k + (j - 1); };
  for (int j = 1; j <= k; ++j) iwork[at(1, j)] = 1;
  for (int i = 1; i <= m; ++i) {
    double ss = 0, s = 0, var = 0;
    for (int ii = 1; ii <= i; ++ii) {
      const int iii = i - ii + 1;
      ss = ss + x[iii - 1] * x[iii - 1];
      s = s + x[iii - 1];
      const double sn = ii;
      var = ss - (s * s) / sn;
      const int ik = iii - 1;
      if (ik != 0) {
        for (int j = 2; j <= k; ++j) {
          const double cand = var + work[at(ik, j - 1)];
          if (work[at(i, j)] >= cand) {
            iwork[at(i, j)] = iii;
            work[at(i, j)] = cand;
          }
        }
      }
    }
    work[at(i, 1)] = var;
    iwork[at(i, 1)] = 1;
  }
  std::vector<int> first(k);
  int il = m + 1;
  for (int l = 1; l <= k; ++l) {
    const int ll = k - l + 1;
    const int iu = il - 1;
    il = iwork[at(iu, ll)];
    first[ll - 1] = il;
  }
  return first;
}

// Same optimum as fisher_exact() in O(k n log n): for a fixed number of classes
// the (smallest) optimal start of the last class is non-decreasing in the
// number of elements, so it can be found by divide and conquer.
class FisherDC {
 public:
  FisherDC(const std::vector<double> &x, const int k) : x_(x), k_(k), n_(static_cast<int>(x.size())) {
    // prefix sums of the centred values (long double for accuracy)
    long double mean = 0;
    for (double v : x_) mean += v;
    mean /= n_;
    p1_.assign(n_ + 1, 0.0L);
    p2_.assign(n_ + 1, 0.0L);
    for (int i = 1; i <= n_; ++i) {
      const long double v = static_cast<long double>(x_[i - 1]) - mean;
      p1_[i] = p1_[i - 1] + v;
      p2_[i] = p2_[i - 1] + v * v;
    }
  }

  std::vector<int> solve() {
    const long double inf = std::numeric_limits<long double>::infinity();
    std::vector<long double> prev(n_ + 1, inf), cur(n_ + 1, inf);
    arg_.assign(static_cast<std::size_t>(k_ + 1) * (n_ + 1), 1);
    for (int i = 1; i <= n_; ++i) prev[i] = cost(1, i);
    for (int j = 2; j <= k_; ++j) {
      std::fill(cur.begin(), cur.end(), inf);
      rec(j, j, n_, j, n_, prev, cur);
      std::swap(prev, cur);
    }
    std::vector<int> first(k_);
    int i = n_;
    for (int j = k_; j >= 2; --j) {
      const int t = arg_[static_cast<std::size_t>(j) * (n_ + 1) + i];
      first[j - 1] = t;
      i = t - 1;
    }
    first[0] = 1;
    return first;
  }

 private:
  // within-class sum of squares of x[t..i] (1-based, inclusive)
  long double cost(const int t, const int i) const {
    const long double s1 = p1_[i] - p1_[t - 1];
    const long double s2 = p2_[i] - p2_[t - 1];
    const long double c = s2 - s1 * s1 / (i - t + 1);
    return c > 0 ? c : 0.0L;
  }

  void rec(const int j, const int lo, const int hi, const int optlo, const int opthi,
           const std::vector<long double> &prev, std::vector<long double> &cur) {
    if (lo > hi) return;
    const int mid = lo + (hi - lo) / 2;
    long double best = std::numeric_limits<long double>::infinity();
    int bt = std::max(optlo, j);
    const int tmax = std::min(mid, opthi);
    for (int t = std::max(optlo, j); t <= tmax; ++t) {
      const long double c = prev[t - 1] + cost(t, mid);
      if (c < best) {  // strict: the smallest start wins ties
        best = c;
        bt = t;
      }
    }
    cur[mid] = best;
    arg_[static_cast<std::size_t>(j) * (n_ + 1) + mid] = bt;
    rec(j, lo, mid - 1, optlo, bt, prev, cur);
    rec(j, mid + 1, hi, bt, opthi, prev, cur);
  }

  const std::vector<double> &x_;
  const int k_, n_;
  std::vector<long double> p1_, p2_;
  std::vector<int> arg_;
};

}  // namespace cgei

#endif
