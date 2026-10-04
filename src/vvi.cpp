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

// What run_vvi() computes:
//   kLists   visible cells and potential viewshed (R cell numbers, sorted) of
//            every observer (plus their sizes)
//   kCounts  only the sizes of both sets per observer
//   kCells   per raster cell the number of observers that see it / whose
//            potential viewshed contains it (`visible_count`, `seen_count`:
//            ncell values each, visible_count zero-initialised; n_seen unset)
// The potential viewshed are the cells reached by a line of sight (static
// geometry) inside the raster with a valid height, plus the observer cell.
enum class VviOut { kLists, kCounts, kCells };

struct VviResult {
  std::vector<std::vector<int>> visible, seen;
  std::vector<int> n_visible, n_seen;
};

VviResult run_vvi(const NumericVector &dsm, const NumericVector &dsm_values, const IntegerVector &x0,
                  const IntegerVector &y0, const NumericVector &h0, const int radius,
                  const int ncores, const bool display_progress, const bool early_stop,
                  const VviOut what, int *visible_count = nullptr, int *seen_count = nullptr) {
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
  const cgei::CoverRuns cov(T);

  const double *dsm_v = dsm_values.begin();
  const double *h0_v = h0.begin();
  const int nthreads = cgei::resolve_threads(ncores);
  // Block maxima for the early-termination bounds (an empty grid if disabled)
  const cgei::BlockMax bm(dsm_v, early_stop ? ras.nrow : 0, early_stop ? ras.ncol : 0, nthreads);
  // valid-cell counts for the size of the potential viewshed (kCounts only)
  const cgei::ValidPrefix valid(dsm_v, what == VviOut::kCounts ? ras.nrow : 0, ras.ncol, nthreads);

  const bool lists = what == VviOut::kLists;
  const bool cells = what == VviOut::kCells;
  VviResult res;
  res.n_visible.assign(n, 0);
  res.n_seen.assign(n, 0);
  if (lists) {
    res.visible.resize(n);
    res.seen.resize(n);
  }

  std::vector<cgei::SweepScratch> scratch;
  scratch.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t) scratch.emplace_back(T);

  const int nvalid = static_cast<int>(obs.order.size());
  cgei::ParallelProgress progress(nvalid, display_progress);

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
    std::vector<int> *vis = lists ? &res.visible[k] : nullptr;
    if (vis) vis->push_back(static_cast<int>(cell0) + 1);
    if (cells) {
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++visible_count[cell0];
    }
    if (hk > dsm_v[cell0]) {
      double qbound[4] = {0, 0, 0, 0};
      if (early_stop) cgei::quadrant_bounds(T, bm, row0, col0, hk, qbound);
      auto visit = [&](const int s, const long long cell) {
        if (!S.first_visit(T.step[s].ref)) return;
        ++n_vis;
        if (vis) vis->push_back(static_cast<int>(cell) + 1);  // C++ -> R cell number
        if (cells) {
#ifdef _OPENMP
#pragma omp atomic
#endif
          ++visible_count[cell];
        }
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

    // Potential viewshed
    if (lists) {
      std::vector<int> &seen = res.seen[k];
      bool own_done = false;
      cov.for_each(row0, col0, ras.nrow, ras.ncol, [&](const int rr, const int c0, const int c1) {
        if (!own_done && (rr > row0 || (rr == row0 && c0 > col0))) {
          // keep the output sorted: the observer cell goes before the first
          // covered cell that follows it in row-major order
          seen.push_back(static_cast<int>(cell0) + 1);
          own_done = true;
        }
        const long long base = static_cast<long long>(rr) * ras.ncol;
        for (int c = c0; c <= c1; ++c) {
          if (!std::isnan(dsm_v[base + c])) seen.push_back(static_cast<int>(base + c) + 1);
        }
      });
      if (!own_done) seen.push_back(static_cast<int>(cell0) + 1);
      res.n_seen[k] = static_cast<int>(seen.size());
    } else if (!cells) {
      res.n_seen[k] = cgei::count_seen(cov, valid, row0, col0, ras.nrow, ras.ncol);
    }
    progress.tick();
  }
  progress.finish();  // throws if the user interrupted

  if (cells) {
    std::vector<int> rows, cols;
    rows.reserve(nvalid);
    cols.reserve(nvalid);
    for (int idx = 0; idx < nvalid; ++idx) {
      rows.push_back(obs.row[obs.order[idx]]);
      cols.push_back(obs.col[obs.order[idx]]);
    }
    cgei::accumulate_seen(cov, dsm_v, ras.nrow, ras.ncol, rows, cols, nthreads, seen_count);
  }
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
                          early_stop, VviOut::kLists);
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
                          early_stop, VviOut::kCounts);
  return Rcpp::List::create(Rcpp::Named("n_visible") = Rcpp::wrap(res.n_visible),
                            Rcpp::Named("n_viewshed") = Rcpp::wrap(res.n_seen));
}

// Per raster cell (R cell order): the number of observers that see the cell
// and the number of observers whose potential viewshed contains it, i.e.
// tabulate() of the cells returned by VVI_cpp(), without the per-observer
// lists (vvi(mode = "cumulative" / "viewshed")).
// [[Rcpp::export]]
Rcpp::List VVI_cells_cpp(const Rcpp::NumericVector &dsm, const Rcpp::NumericVector &dsm_values,
                         const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                         const Rcpp::NumericVector &h0, const int radius,
                         const int ncores = 1, const bool display_progress = false,
                         const bool early_stop = true) {
  const RasterInfo ras(dsm);
  const R_xlen_t ncell = static_cast<R_xlen_t>(ras.nrow) * ras.ncol;
  if (static_cast<double>(ncell) >= 2147483647.0)
    Rcpp::stop("The DSM has too many cells to return R cell numbers.");
  Rcpp::IntegerVector visible_count(ncell);  // zero-initialised
  Rcpp::IntegerVector seen_count(ncell);
  // raw pointers: no R object is touched inside the parallel region
  run_vvi(dsm, dsm_values, x0, y0, h0, radius, ncores, display_progress, early_stop,
          VviOut::kCells, visible_count.begin(), seen_count.begin());
  return Rcpp::List::create(Rcpp::Named("visible_count") = visible_count,
                            Rcpp::Named("viewshed_count") = seen_count);
}
