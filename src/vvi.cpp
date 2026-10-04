#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "parallel_progress.h"
#include "rsinfo.h"
#include "viewshed_engine.h"

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

struct VviObservers {
  std::vector<long long> cell;  // 0-based linear index, -1 if invalid
  std::vector<int> row, col;
  std::vector<int> order;       // valid observers in Morton order
};

VviObservers prepare(const IntegerVector &x0, const IntegerVector &y0, const RasterInfo &ras) {
  if (x0.size() != y0.size()) Rcpp::stop("x0 and y0 must have the same length.");
  const int n = x0.size();
  VviObservers o;
  o.cell.assign(n, -1);
  o.row.assign(n, -1);
  o.col.assign(n, -1);
  std::vector<std::pair<std::uint64_t, int>> keys;
  keys.reserve(n);
  for (int k = 0; k < n; ++k) {
    if (IntegerVector::is_na(x0[k]) || IntegerVector::is_na(y0[k])) continue;
    const int col = x0[k] - 1;  // R (1-based) -> C++ (0-based)
    const int row = y0[k] - 1;
    if (row < 0 || row >= ras.nrow || col < 0 || col >= ras.ncol) continue;
    o.row[k] = row;
    o.col[k] = col;
    o.cell[k] = static_cast<long long>(row) * ras.ncol + col;
    keys.emplace_back(cgei::morton_key(row, col), k);
  }
  std::sort(keys.begin(), keys.end());
  for (const auto &p : keys) o.order.push_back(p.second);
  return o;
}

// Visible cells (1-based, sorted) and potential viewshed (cells reached by a
// line of sight inside the raster with a valid height; 1-based, sorted) of
// every observer. If counts_only, only the two set sizes are stored.
struct VviResult {
  std::vector<std::vector<int>> visible, seen;
  std::vector<int> n_visible, n_seen;
};

VviResult run_vvi(const NumericVector &dsm, const NumericVector &dsm_values, const IntegerVector &x0,
                  const IntegerVector &y0, const NumericVector &h0, const int radius,
                  const int ncores, const bool display_progress, const bool early_stop,
                  const bool counts_only) {
  const RasterInfo ras(dsm);
  if (dsm_values.size() != static_cast<R_xlen_t>(ras.nrow) * ras.ncol)
    Rcpp::stop("dsm_values does not match the dimensions of dsm.");
  if (static_cast<double>(ras.nrow) * ras.ncol >= 2147483647.0)
    Rcpp::stop("The DSM has too many cells to return R cell numbers.");
  if (radius < 1) Rcpp::stop("radius must be at least 1.");

  const VviObservers obs = prepare(x0, y0, ras);
  const int n = static_cast<int>(obs.cell.size());
  if (h0.size() != n) Rcpp::stop("h0 must have the same length as x0.");

  const int r = static_cast<int>(std::round(radius / ras.res));
  cgei::LosTable T(r);
  if (!T.bind(ras.ncol)) Rcpp::stop("Raster too large for the radius.");

  const double *dsm_v = dsm_values.begin();
  const double *h0_v = h0.begin();
  const int nthreads = cgei::resolve_threads(ncores);
  // Block maxima for the early-termination bounds (an empty grid if disabled)
  const cgei::BlockMax bm(dsm_v, early_stop ? ras.nrow : 0, early_stop ? ras.ncol : 0, nthreads);

  VviResult res;
  res.n_visible.assign(n, 0);
  res.n_seen.assign(n, 0);
  if (!counts_only) {
    res.visible.resize(n);
    res.seen.resize(n);
  }

  std::vector<cgei::SweepScratch> scratch;
  scratch.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t) scratch.emplace_back(T);

  const int nvalid = static_cast<int>(obs.order.size());
  cgei::ParallelProgress progress(nvalid, display_progress);
  const int ncov = static_cast<int>(T.cover_dr.size());

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
#endif
  for (int idx = 0; idx < nvalid; ++idx) {
    if (progress.aborted()) continue;
    cgei::SweepScratch &S = scratch[cgei::omp_thread_num()];
    const int k = obs.order[idx];
    const long long cell0 = obs.cell[k];
    const int row0 = obs.row[k], col0 = obs.col[k];
    const double hk = h0_v[k];

    // Visible cells: the observer cell plus everything found by the sweep.
    S.first_visit(T.ref_center());
    int n_vis = 1;
    std::vector<int> *vis = counts_only ? nullptr : &res.visible[k];
    if (vis) vis->push_back(static_cast<int>(cell0) + 1);
    if (hk > dsm_v[cell0]) {
      double qbound[4] = {0, 0, 0, 0};
      if (early_stop) cgei::quadrant_bounds(T, bm, row0, col0, hk, qbound);
      auto visit = [&](const int s, const long long cell) {
        if (!S.first_visit(T.step[s].ref)) return;
        ++n_vis;
        if (vis) vis->push_back(static_cast<int>(cell) + 1);  // C++ -> R cell number
      };
      const bool interior = row0 - r >= 0 && row0 + r < ras.nrow && col0 - r >= 0 &&
                            col0 + r < ras.ncol;
      if (interior) {
        cgei::sweep_lines<false>(T, dsm_v, ras.nrow, ras.ncol, cell0, row0, col0, hk,
                                 qbound, early_stop, S, visit);
      } else {
        cgei::sweep_lines<true>(T, dsm_v, ras.nrow, ras.ncol, cell0, row0, col0, hk,
                                qbound, early_stop, S, visit);
      }
    }
    S.reset();
    if (vis) std::sort(vis->begin(), vis->end());
    res.n_visible[k] = n_vis;

    // Potential viewshed: all cells reached by any line of sight (static
    // geometry) that lie inside the raster and have a valid height.
    std::vector<int> *seen = counts_only ? nullptr : &res.seen[k];
    int n_seen = 0;
    bool own_done = false;
    for (int c = 0; c < ncov; ++c) {
      const int rr = row0 + T.cover_dr[c];
      const int cc = col0 + T.cover_dc[c];
      if (!own_done && (rr > row0 || (rr == row0 && cc > col0))) {
        // keep the output sorted: the observer cell goes before the first
        // covered cell that follows it in row-major order
        ++n_seen;
        if (seen) seen->push_back(static_cast<int>(cell0) + 1);
        own_done = true;
      }
      if (rr < 0 || rr >= ras.nrow || cc < 0 || cc >= ras.ncol) continue;
      const long long cell = static_cast<long long>(rr) * ras.ncol + cc;
      if (std::isnan(dsm_v[cell])) continue;
      ++n_seen;
      if (seen) seen->push_back(static_cast<int>(cell) + 1);
    }
    if (!own_done) {
      ++n_seen;
      if (seen) seen->push_back(static_cast<int>(cell0) + 1);
    }
    res.n_seen[k] = n_seen;
    progress.tick();
  }
  progress.finish();  // throws if the user interrupted
  return res;
}

}  // namespace

// Visible cells and potential viewshed of every observer (R cell numbers).
// [[Rcpp::export]]
Rcpp::List VVI_cpp(const Rcpp::NumericVector &dsm, const Rcpp::NumericVector &dsm_values,
                   const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                   const Rcpp::NumericVector &h0, const int radius,
                   const int ncores = 1, const bool display_progress = false,
                   const bool early_stop = true) {
  VviResult res = run_vvi(dsm, dsm_values, x0, y0, h0, radius, ncores, display_progress,
                          early_stop, false);
  const int n = static_cast<int>(res.n_visible.size());
  Rcpp::List output(n);
  for (int k = 0; k < n; ++k) {
    output[k] = Rcpp::List::create(Rcpp::Named("visible_cells") = Rcpp::wrap(res.visible[k]),
                                   Rcpp::Named("viewshed") = Rcpp::wrap(res.seen[k]));
  }
  return output;
}

// Only the number of visible cells and of potential viewshed cells per
// observer (enough for vvi(mode = "VVI")); 0 for invalid observers.
// [[Rcpp::export]]
Rcpp::List VVI_count_cpp(const Rcpp::NumericVector &dsm, const Rcpp::NumericVector &dsm_values,
                         const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                         const Rcpp::NumericVector &h0, const int radius,
                         const int ncores = 1, const bool display_progress = false,
                         const bool early_stop = true) {
  VviResult res = run_vvi(dsm, dsm_values, x0, y0, h0, radius, ncores, display_progress,
                          early_stop, true);
  return Rcpp::List::create(Rcpp::Named("n_visible") = Rcpp::wrap(res.n_visible),
                            Rcpp::Named("n_viewshed") = Rcpp::wrap(res.n_seen));
}
