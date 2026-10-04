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

// Everything the parallel loop needs (plain C++ objects and raw pointers only).
struct VviJob {
  const RasterInfo *ras;
  const VviObservers *obs;
  const cgei::LosTable *T;
  const cgei::VisMask *M;
  const cgei::CoverRuns *cov;
  const cgei::BlockMax *bm;
  const cgei::ValidPrefix *valid;
  const double *dsm_v;
  const double *h0_v;
  int r, nthreads, batch;
  bool early_stop;
  VviOut what;
  int *visible_count;
};

// The observers in Morton order, in batches of job.batch: sweep a batch, then
// read every observer's mask (in raster order, i.e. sorted by cell number).
// V is the storage type of the DSM used by the sweep (float if exact).
template <class V>
void vvi_loop(const VviJob &job, const V *dsm_sweep, cgei::ParallelProgress &progress,
              VviResult &res) {
  const RasterInfo &ras = *job.ras;
  const VviObservers &obs = *job.obs;
  const cgei::LosTable &T = *job.T;
  const cgei::VisMask &M = *job.M;
  const int r = job.r, nrow = ras.nrow, ncol = ras.ncol;
  const bool lists = job.what == VviOut::kLists;
  const bool cells = job.what == VviOut::kCells;

  std::vector<cgei::BatchScratch> scratch;
  scratch.reserve(job.nthreads);
  for (int t = 0; t < job.nthreads; ++t) scratch.emplace_back(T, M, job.batch);

  const int nvalid = static_cast<int>(obs.order.size());
  const int nbatch = (nvalid + job.batch - 1) / job.batch;

#ifdef _OPENMP
#pragma omp parallel for num_threads(job.nthreads) schedule(dynamic, 1)
#endif
  for (int b = 0; b < nbatch; ++b) {
    if (progress.aborted()) continue;
    const int tid = cgei::omp_thread_num();
    cgei::BatchScratch &S = scratch[tid];
    const int i0 = b * job.batch, i1 = std::min(nvalid, i0 + job.batch);

    // Visible cells: the observer cell plus everything found by the sweep.
    S.viewers.clear();
    for (int i = i0; i < i1; ++i) {
      const int k = obs.order[i];
      std::uint64_t *mask = S.mask(i - i0);
      cgei::mask_set(mask, M.center);
      if (job.h0_v[k] > job.dsm_v[obs.cell[k]]) {
        cgei::Viewer v;
        v.cell0 = obs.cell[k];
        v.row0 = obs.row[k];
        v.col0 = obs.col[k];
        v.h0 = job.h0_v[k];
        v.interior = v.row0 - r >= 0 && v.row0 + r < nrow && v.col0 - r >= 0 && v.col0 + r < ncol;
        v.mask = mask;
        v.hz = S.horizon(i - i0);
        v.valid_upto = -1;
        if (job.early_stop) cgei::quadrant_bounds(T, *job.bm, v.row0, v.col0, v.h0, v.qbound);
        S.viewers.push_back(v);
      }
    }
    cgei::sweep_batch(T, M, dsm_sweep, nrow, ncol, job.early_stop, S.viewers.data(),
                      static_cast<int>(S.viewers.size()));

    for (int i = i0; i < i1; ++i) {
      const int k = obs.order[i];
      const int row0 = obs.row[k], col0 = obs.col[k];
      const long long cell0 = obs.cell[k];
      std::uint64_t *mask = S.mask(i - i0);
      if (lists) {
        std::vector<int> &vis = res.visible[k];
        cgei::drain_mask(M, mask, [&](const int dr, const int dc) {
          vis.push_back(static_cast<int>(cell0 + static_cast<long long>(dr) * ncol + dc) + 1);  // C++ -> R
        });
        res.n_visible[k] = static_cast<int>(vis.size());

        // Potential viewshed (sorted: runs are row-major, the observer cell is
        // inserted before the first covered cell that follows it)
        std::vector<int> &seen = res.seen[k];
        bool own_done = false;
        job.cov->for_each(row0, col0, nrow, ncol, [&](const int rr, const int c0, const int c1) {
          if (!own_done && (rr > row0 || (rr == row0 && c0 > col0))) {
            seen.push_back(static_cast<int>(cell0) + 1);
            own_done = true;
          }
          const long long base = static_cast<long long>(rr) * ncol;
          for (int c = c0; c <= c1; ++c) {
            if (!std::isnan(job.dsm_v[base + c])) seen.push_back(static_cast<int>(base + c) + 1);
          }
        });
        if (!own_done) seen.push_back(static_cast<int>(cell0) + 1);
        res.n_seen[k] = static_cast<int>(seen.size());
      } else if (cells) {
        int n_vis = 0;
        cgei::drain_mask(M, mask, [&](const int dr, const int dc) {
          const long long cell = cell0 + static_cast<long long>(dr) * ncol + dc;
#ifdef _OPENMP
#pragma omp atomic
#endif
          ++job.visible_count[cell];
          ++n_vis;
        });
        res.n_visible[k] = n_vis;
      } else {
        res.n_visible[k] = cgei::count_mask(M, mask);
        res.n_seen[k] = cgei::count_seen(*job.cov, *job.valid, row0, col0, nrow, ncol);
      }
    }
    progress.tick(static_cast<unsigned long>(i1 - i0));
  }
}

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
  if (r > cgei::kMaxRadius) Rcpp::stop("max_distance is too large for the resolution of the DSM.");
  cgei::LosTable T(r);
  if (!T.bind(ras.ncol)) Rcpp::stop("Raster too large for the radius.");
  const cgei::VisMask M(r);
  const cgei::CoverRuns cov(T);

  const double *dsm_v = dsm_values.begin();
  const int nthreads = cgei::resolve_threads(ncores);
  // Block maxima for the early-termination bounds (an empty grid if disabled)
  const cgei::BlockMax bm(dsm_v, early_stop ? ras.nrow : 0, early_stop ? ras.ncol : 0, nthreads);
  // valid-cell counts for the size of the potential viewshed (kCounts only)
  const cgei::ValidPrefix valid(dsm_v, what == VviOut::kCounts ? ras.nrow : 0, ras.ncol, nthreads);

  VviResult res;
  res.n_visible.assign(n, 0);
  res.n_seen.assign(n, 0);
  if (what == VviOut::kLists) {
    res.visible.resize(n);
    res.seen.resize(n);
  }

  VviJob job;
  job.ras = &ras;
  job.obs = &obs;
  job.T = &T;
  job.M = &M;
  job.cov = &cov;
  job.bm = &bm;
  job.valid = &valid;
  job.dsm_v = dsm_v;
  job.h0_v = h0.begin();
  job.r = r;
  job.nthreads = nthreads;
  job.batch = cgei::batch_size(M, static_cast<int>(obs.order.size()), nthreads);
  job.early_stop = early_stop;
  job.what = what;
  job.visible_count = visible_count;

  cgei::ParallelProgress progress(obs.order.size(), display_progress);
  const std::size_t ncell = static_cast<std::size_t>(ras.nrow) * ras.ncol;
  if (cgei::float_exact(dsm_v, ncell, nthreads)) {
    std::vector<float> dsm_f(ncell);
    cgei::to_float(dsm_v, ncell, dsm_f.data(), nthreads);
    vvi_loop<float>(job, dsm_f.data(), progress, res);
  } else {
    vvi_loop<double>(job, dsm_v, progress, res);
  }
  progress.finish();  // throws if the user interrupted

  if (what == VviOut::kCells) {
    const int nvalid = static_cast<int>(obs.order.size());
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
