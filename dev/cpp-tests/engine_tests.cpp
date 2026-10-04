// Stand-alone tests of CGEI's C++ engines (no R needed), meant to be run under
// sanitizers (see Makefile):
//   make asan   AddressSanitizer + UndefinedBehaviorSanitizer (gcc)
//   make tsan   ThreadSanitizer with the LLVM OpenMP runtime + Archer (clang)
//
// The headers tested here are exactly the ones compiled into the R package.
// Every engine is run multi-threaded and compared with a naive, single-threaded
// re-implementation written independently below.

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#include "../../src/boxfilter.h"
#include "../../src/los_geometry.h"
#include "../../src/natural_breaks.h"
#include "../../src/viewshed_engine.h"

using namespace cgei;

static int n_checks = 0, n_failures = 0;
#define CHECK(cond, ...)                                                    \
  do {                                                                      \
    ++n_checks;                                                             \
    if (!(cond)) {                                                          \
      ++n_failures;                                                         \
      std::fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__);             \
      std::fprintf(stderr, __VA_ARGS__);                                    \
      std::fprintf(stderr, "\n");                                           \
    }                                                                       \
  } while (0)

static const double NaN = std::numeric_limits<double>::quiet_NaN();

std::vector<double> random_dsm(int nr, int nc, double na_frac, std::mt19937 &rng) {
  std::uniform_real_distribution<double> u(0, 1);
  std::vector<double> m(static_cast<size_t>(nr) * nc);
  for (auto &v : m) v = 0.5 * u(rng);
  std::uniform_int_distribution<int> rr(0, nr - 1), cc(0, nc - 1), sz(0, 3);
  for (int b = 0; b < nr * nc / 60 + 1; ++b) {  // buildings / trees
    const int r0 = rr(rng), c0 = cc(rng), h = 2 + static_cast<int>(18 * u(rng));
    const int r1 = std::min(nr - 1, r0 + sz(rng)), c1 = std::min(nc - 1, c0 + sz(rng));
    for (int i = r0; i <= r1; ++i)
      for (int j = c0; j <= c1; ++j) m[static_cast<size_t>(i) * nc + j] = h;
  }
  for (auto &v : m)
    if (u(rng) < na_frac) v = NaN;
  return m;
}

// ---------------------------------------------------------------------------
// 1. Viewshed engine vs naive line walk
// ---------------------------------------------------------------------------
std::vector<long long> naive_visible(const std::vector<int> &los, int r, const double *dsm,
                                     int nr, int nc, int row0, int col0, double h0) {
  std::vector<long long> vis;
  const long long cell0 = static_cast<long long>(row0) * nc + col0;
  vis.push_back(cell0);
  if (h0 > dsm[cell0]) {
    const int ncr = 2 * r + 1;
    for (int l = 0; l < 8 * r; ++l) {
      double horizon = -9999.0;
      for (int j = 0; j < r; ++j) {
        const int ref = los[static_cast<size_t>(l) * r + j];
        if (ref == kNaInt) break;
        const int dr = ref / ncr - r, dc = ref % ncr - r;
        const int rr = row0 + dr, cc = col0 + dc;
        if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
        const double h = dsm[static_cast<size_t>(rr) * nc + cc];
        if (std::isnan(h)) continue;
        const double t = (h - h0) / std::sqrt(static_cast<double>(dr * dr + dc * dc));
        if (t > horizon) {
          horizon = t;
          vis.push_back(static_cast<long long>(rr) * nc + cc);
        }
      }
    }
  }
  std::sort(vis.begin(), vis.end());
  vis.erase(std::unique(vis.begin(), vis.end()), vis.end());
  return vis;
}

void test_viewshed(std::mt19937 &rng, int nthreads) {
  struct Scene { int nr, nc, r; double na; };
  const Scene scenes[] = {{40, 55, 9, 0.0}, {40, 55, 9, 0.05}, {35, 11, 9, 0.02},
                          {1, 30, 5, 0.0}, {30, 1, 5, 0.0}, {25, 25, 1, 0.0},
                          {60, 45, 17, 0.03}, {12, 12, 2, 0.1}};
  for (const Scene &sc : scenes) {
    std::vector<double> dsm = random_dsm(sc.nr, sc.nc, sc.na, rng);
    LosTable T(sc.r);
    CHECK(T.bind(sc.nc), "bind");
    const std::vector<int> los = los_reference(sc.r, sc.r, sc.r, 2 * sc.r + 1);
    BlockMax bm(dsm.data(), sc.nr, sc.nc, nthreads);
    // observers: every 3rd cell incl. all borders
    std::vector<int> rows, cols;
    std::vector<double> h0;
    std::uniform_real_distribution<double> u(0, 3);
    for (int i = 0; i < sc.nr; ++i)
      for (int j = 0; j < sc.nc; ++j)
        if ((i * sc.nc + j) % 3 == 0 || i == 0 || j == 0 || i == sc.nr - 1 || j == sc.nc - 1) {
          rows.push_back(i);
          cols.push_back(j);
          h0.push_back(u(rng));
        }
    const int n = static_cast<int>(rows.size());
    for (int early = 0; early <= 1; ++early) {
      std::vector<std::vector<long long>> out(n);
      std::vector<SweepScratch> scratch;
      for (int t = 0; t < nthreads; ++t) scratch.emplace_back(T);
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 4)
      for (int k = 0; k < n; ++k) {
        SweepScratch &S = scratch[omp_get_thread_num()];
        const long long cell0 = static_cast<long long>(rows[k]) * sc.nc + cols[k];
        out[k].push_back(cell0);
        S.first_visit(T.ref_center());
        if (h0[k] > dsm[cell0]) {
          double qb[4] = {0, 0, 0, 0};
          if (early) quadrant_bounds(T, bm, rows[k], cols[k], h0[k], qb);
          auto visit = [&](int s, long long cell) {
            if (S.first_visit(T.step[s].ref)) out[k].push_back(cell);
          };
          sweep_lines<true>(T, dsm.data(), sc.nr, sc.nc, cell0, rows[k], cols[k], h0[k], qb,
                            early == 1, S, visit);
        }
        S.reset();
        std::sort(out[k].begin(), out[k].end());
      }
      int bad = 0;
      for (int k = 0; k < n; ++k) {
        if (std::isnan(dsm[static_cast<size_t>(rows[k]) * sc.nc + cols[k]])) continue;
        if (out[k] != naive_visible(los, sc.r, dsm.data(), sc.nr, sc.nc, rows[k], cols[k], h0[k])) ++bad;
      }
      CHECK(bad == 0, "viewshed %dx%d r=%d na=%.2f early=%d threads=%d: %d of %d observers differ",
            sc.nr, sc.nc, sc.r, sc.na, early, nthreads, bad, n);
    }
  }
}

// ---------------------------------------------------------------------------
// 2. Box filters vs naive loops
// ---------------------------------------------------------------------------
void test_boxfilters(std::mt19937 &rng, int nthreads) {
  std::uniform_int_distribution<int> val(0, 9);
  std::uniform_real_distribution<double> u(0, 1);
  const int shapes[][2] = {{17, 23}, {1, 40}, {40, 1}, {5, 5}, {64, 129}};
  for (const auto &sh : shapes) {
    const int nr = sh[0], nc = sh[1];
    std::vector<double> v(static_cast<size_t>(nr) * nc);
    for (auto &x : v) x = (u(rng) < 0.1) ? NaN : val(rng);  // integer data: exact sums
    // centred windows (gavi): sums and counts, clipped at the borders
    for (int rad : {0, 1, 3, 10, 70}) {
      std::vector<double> H(v.size()), S(v.size());
      std::vector<int> Hc(v.size()), C(v.size());
      box_rows(v.data(), nr, nc, rad, rad, nc, H.data(), Hc.data(), nthreads);
      box_cols(H.data(), Hc.data(), nr, nc, rad, rad, nr, nthreads,
               [&](int row, int c0, int c1, const double *sum, const int *cnt) {
                 for (int c = c0; c < c1; ++c) {
                   S[static_cast<size_t>(row) * nc + c] = sum[c - c0];
                   C[static_cast<size_t>(row) * nc + c] = cnt[c - c0];
                 }
               });
      int bad = 0;
      for (int i = 0; i < nr; ++i)
        for (int j = 0; j < nc; ++j) {
          double s = 0;
          int n = 0;
          for (int a = std::max(0, i - rad); a <= std::min(nr - 1, i + rad); ++a)
            for (int b = std::max(0, j - rad); b <= std::min(nc - 1, j + rad); ++b) {
              const double x = v[static_cast<size_t>(a) * nc + b];
              if (!std::isnan(x)) {
                s += x;
                ++n;
              }
            }
          if (S[static_cast<size_t>(i) * nc + j] != s || C[static_cast<size_t>(i) * nc + j] != n) ++bad;
        }
      CHECK(bad == 0, "centred box sums %dx%d rad=%d threads=%d: %d cells differ", nr, nc, rad, nthreads, bad);
    }
    // anchored w x w boxes (lacunarity): sums, max, min
    for (int w : {1, 2, 3, 5, 16}) {
      if (w > nr || w > nc) continue;
      const int onr = nr - w + 1, onc = nc - w + 1;
      std::vector<double> H(static_cast<size_t>(nr) * onc), S(static_cast<size_t>(onr) * onc);
      std::vector<double> mx(S.size()), mn(S.size());
      box_rows(v.data(), nr, nc, 0, w - 1, onc, H.data(), nullptr, nthreads);
      box_cols(H.data(), nullptr, nr, onc, 0, w - 1, onr, nthreads,
               [&](int row, int c0, int c1, const double *sum, const int *) {
                 for (int c = c0; c < c1; ++c) S[static_cast<size_t>(row) * onc + c] = sum[c - c0];
               });
      box_extreme<true>(v.data(), nr, nc, w, mx.data(), nthreads);
      box_extreme<false>(v.data(), nr, nc, w, mn.data(), nthreads);
      int bad = 0;
      for (int i = 0; i < onr; ++i)
        for (int j = 0; j < onc; ++j) {
          double s = 0, hi = -std::numeric_limits<double>::infinity(),
                 lo = std::numeric_limits<double>::infinity();
          for (int a = i; a < i + w; ++a)
            for (int b = j; b < j + w; ++b) {
              const double x = v[static_cast<size_t>(a) * nc + b];
              if (std::isnan(x)) continue;
              s += x;
              hi = std::max(hi, x);
              lo = std::min(lo, x);
            }
          const size_t o = static_cast<size_t>(i) * onc + j;
          if (S[o] != s || mx[o] != hi || mn[o] != lo) ++bad;
        }
      CHECK(bad == 0, "anchored boxes %dx%d w=%d threads=%d: %d boxes differ", nr, nc, w, nthreads, bad);
    }
  }
}

// ---------------------------------------------------------------------------
// 3. Natural breaks: O(k n log n) Fisher == exact O(k n^2) Fortran port
// ---------------------------------------------------------------------------
// within-class sum of squares of a partition (long double, from scratch)
long double partition_sse(const std::vector<double> &x, const std::vector<int> &first) {
  long double total = 0;
  const int k = static_cast<int>(first.size());
  for (int c = 0; c < k; ++c) {
    const int a = first[c] - 1, b = c + 1 < k ? first[c + 1] - 1 : static_cast<int>(x.size());
    long double m = 0;
    for (int i = a; i < b; ++i) m += x[i];
    m /= (b - a);
    for (int i = a; i < b; ++i) total += (x[i] - m) * (x[i] - m);
  }
  return total;
}

void test_natural_breaks(std::mt19937 &rng) {
  std::uniform_real_distribution<double> u(0, 1);
  std::normal_distribution<double> g(0, 1);
  std::uniform_int_distribution<int> d12(0, 12);
  int identical = 0, worse = 0, total = 0, identical_cont = 0, total_cont = 0;
  for (int rep = 0; rep < 60; ++rep) {
    const int n = 20 + rep * 23;
    const int family = rep % 3;  // 0: uniform, 1: normal mixture, 2: 13 distinct values
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
      x[i] = family == 0 ? u(rng) : family == 1 ? (i % 2 ? g(rng) : 4 + g(rng)) : d12(rng) / 12.0;
    }
    std::sort(x.begin(), x.end());
    for (int k = 2; k <= 9; ++k) {
      FisherDC dc(x, k);
      const std::vector<int> a = dc.solve(), b = fisher_exact(x, k);
      ++total;
      if (a == b) ++identical;
      // the fast version may only differ in (near) ties and is never worse
      if (partition_sse(x, a) > partition_sse(x, b) * (1 + 1e-12L) + 1e-18L) ++worse;
      if (family != 2) {
        ++total_cont;
        if (a == b) ++identical_cont;
      }
    }
  }
  std::printf("natural breaks: Fisher D&C identical to the exact port in %d / %d cases "
              "(continuous data: %d / %d), never worse: %s\n",
              identical, total, identical_cont, total_cont, worse == 0 ? "yes" : "NO");
  CHECK(worse == 0, "Fisher D&C found a worse partition than the exact port in %d cases", worse);
  CHECK(identical_cont == total_cont, "Fisher D&C differs on continuous data (%d of %d identical)",
        identical_cont, total_cont);
  std::vector<double> x = {1, 2, 3, 10, 11, 12, 20, 21, 22};
  CHECK((jenks_classint(x, 3) == std::vector<double>{1, 3, 12, 22}), "jenks hand example");
  FisherDC dc(x, 3);
  CHECK((fisher_breaks(x, dc.solve()) == std::vector<double>{1, 6.5, 16, 22}), "fisher hand example");
}

// ---------------------------------------------------------------------------
// 4. Structural port of the OLD VGVI/VVI viewshed loop (CGEI 0.3.1, including
//    its bugs; std::vector instead of Rcpp vectors). Only used to check the
//    OpenMP structure of the original code with ThreadSanitizer.
// ---------------------------------------------------------------------------
std::vector<int> old_loop(const std::vector<double> &dsm, int nrow, int ncol,
                          const std::vector<int> &x0, const std::vector<int> &y0,
                          const std::vector<double> &h0, int r, int nthreads) {
  const int nc_ref = 2 * r + 1, c0_ref = r * nc_ref + r;
  const std::vector<int> los_ref_vec = los_reference(r, r, r, nc_ref);
  const std::vector<int> los_start = shared_los(r, los_ref_vec);
  const int n = static_cast<int>(x0.size());
  std::vector<int> output(n, -1);
#pragma omp parallel for num_threads(nthreads) schedule(dynamic) shared(output)
  for (int k = 0; k < n; ++k) {
    const int this_input_cell = y0[k] * ncol + x0[k];
    std::vector<int> viewshed(static_cast<size_t>(nc_ref) * nc_ref, kNaInt);
    viewshed[c0_ref] = this_input_cell + 1;
    if (h0[k] > dsm[this_input_cell]) {
      const int x = this_input_cell - c0_ref - r * (ncol - nc_ref);
      std::vector<double> max_tan_vec(r, -9999.0);
      for (int i = 0; i < r * 8; ++i) {
        const int k_i = los_start[i];
        double max_tan = (k_i > 1) ? max_tan_vec[k_i - 1] : -9999.0;
        for (int j = k_i; j < r; ++j) {
          const int los_ref_cell = los_ref_vec[static_cast<size_t>(i) * r + j];
          if (los_ref_cell == kNaInt) break;
          const int cell = x + los_ref_cell + (los_ref_cell / nc_ref) * (ncol - nc_ref);
          const int row = cell / ncol, col = cell - row * ncol;
          if (!(cell < 0 || cell >= nrow * ncol || std::abs(col - x0[k]) > r)) {
            const double h_cell = dsm[cell];
            if (std::isnan(h_cell)) continue;
            const double d = std::sqrt(static_cast<double>((x0[k] - col) * (x0[k] - col) +
                                                           (y0[k] - row) * (y0[k] - row)));
            const double this_tan = (h_cell - h0[k]) / d;
            if (this_tan > max_tan) {
              max_tan = this_tan;
              viewshed[los_ref_cell] = cell + 1;
            }
          }
          max_tan_vec[j] = max_tan;
        }
      }
    }
    viewshed.erase(std::remove(viewshed.begin(), viewshed.end(), kNaInt), viewshed.end());
    output[k] = static_cast<int>(viewshed.size());
  }
  return output;
}

void test_old_loop(std::mt19937 &rng) {
  const int nr = 60, nc = 70, r = 10;
  std::vector<double> dsm = random_dsm(nr, nc, 0.03, rng);
  std::vector<int> x0, y0;
  std::vector<double> h0;
  for (int i = r; i < nr - r; i += 2)
    for (int j = r; j < nc - r; j += 2)
      if (!std::isnan(dsm[static_cast<size_t>(i) * nc + j])) {
        y0.push_back(i);
        x0.push_back(j);
        h0.push_back(1.7);
      }
  const std::vector<int> a = old_loop(dsm, nr, nc, x0, y0, h0, r, 1);
  const std::vector<int> b = old_loop(dsm, nr, nc, x0, y0, h0, r, 4);
  CHECK(a == b, "old loop: results differ between 1 and 4 threads");
}

int main() {
  std::mt19937 rng(20241004);
  for (int nthreads : {1, 4}) {
    test_viewshed(rng, nthreads);
    test_boxfilters(rng, nthreads);
  }
  test_natural_breaks(rng);
  test_old_loop(rng);
  std::printf("%d checks, %d failures\n", n_checks, n_failures);
  return n_failures == 0 ? 0 : 1;
}
