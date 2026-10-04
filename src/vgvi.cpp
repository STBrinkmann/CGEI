#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "integrate.h"
#include "parallel_progress.h"
#include "rsinfo.h"
#include "viewshed_engine.h"

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

// Observer cells. x0 / y0 are R's 1-based column / row numbers
// (terra::colFromX() / terra::rowFromY()); NA or outside the raster -> invalid.
struct Observers {
  std::vector<long long> cell;  // 0-based linear cell index, -1 if invalid
  std::vector<int> row, col;    // 0-based
  std::vector<int> order;       // valid observers, in Morton (Z-curve) order
};

Observers prepare_observers(const IntegerVector &x0, const IntegerVector &y0,
                            const RasterInfo &ras) {
  if (x0.size() != y0.size()) Rcpp::stop("x0 and y0 must have the same length.");
  const int n = x0.size();
  Observers o;
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
  o.order.reserve(keys.size());
  for (const auto &p : keys) o.order.push_back(p.second);
  return o;
}

// Maps DSM rows / columns to greenspace rows / columns. Replicates the
// original cellFromXY() evaluated at the DSM cell centres, i.e. it also works
// if the greenspace raster has a different extent or resolution.
struct GreenspaceMap {
  std::vector<int> row_of, col_of;
  int gs_ncol = 0;
  bool same_grid = false;  // greenspace grid identical to the DSM grid

  GreenspaceMap(const RasterInfo &dsm, const RasterInfo &gs) : gs_ncol(gs.ncol) {
    same_grid = gs.nrow == dsm.nrow && gs.ncol == dsm.ncol && gs.xmin == dsm.xmin &&
                gs.xmax == dsm.xmax && gs.ymin == dsm.ymin && gs.ymax == dsm.ymax;
    const double yr_inv = gs.nrow / (gs.ymax - gs.ymin);
    const double xr_inv = gs.ncol / (gs.xmax - gs.xmin);
    row_of.resize(dsm.nrow);
    for (int row = 0; row < dsm.nrow; ++row) {
      const double y = dsm.ymax - (row + 0.5) * dsm.res;
      double gr = std::floor((gs.ymax - y) * yr_inv);
      if (y == gs.ymin) gr = gs.nrow - 1;
      row_of[row] = (gr >= 0 && gr < gs.nrow) ? static_cast<int>(gr) : -1;
    }
    col_of.resize(dsm.ncol);
    for (int col = 0; col < dsm.ncol; ++col) {
      const double x = dsm.xmin + (col + 0.5) * dsm.res;
      double gc = std::floor((x - gs.xmin) * xr_inv);
      if (x == gs.xmax) gc = gs.ncol - 1;
      col_of[col] = (gc >= 0 && gc < gs.ncol) ? static_cast<int>(gc) : -1;
    }
  }

  // Greenspace value at the DSM cell (row, col) with linear index cell;
  // NA or outside the greenspace raster -> 0.
  inline double value(const double *gsv, const int row, const int col,
                      const long long cell) const {
    double g;
    if (same_grid) {
      g = gsv[cell];
    } else {
      const int gr = row_of[row], gc = col_of[col];
      if (gr < 0 || gc < 0) return 0.0;
      g = gsv[static_cast<std::size_t>(gr) * gs_ncol + gc];
    }
    return std::isnan(g) ? 0.0 : g;
  }
};

// Per-thread distance-ring histogram (padded, see cgei::kCachePad).
struct Rings {
  char pad_front_[cgei::kCachePad];
  std::vector<int> total;      // visible cells per ring
  std::vector<double> green;   // summed greenspace values per ring
  int max_used = 0;
  char pad_back_[cgei::kCachePad];
  explicit Rings(const int max_ring) : total(max_ring + 1, 0), green(max_ring + 1, 0.0) {}
  inline void add(const int ring, const double g) {
    total[ring] += 1;
    green[ring] += g;
    if (ring > max_used) max_used = ring;
  }
  inline void reset() {
    for (int i = 0; i <= max_used; ++i) {
      total[i] = 0;
      green[i] = 0.0;
    }
    max_used = 0;
  }
};

// VGVI of one observer: (decay weighted) mean over the distance rings that
// contain at least one visible cell of the proportion of green visible cells.
double vgvi_index(const Rings &R, const int fun, const std::vector<double> &weight) {
  if (fun == 3) {  // no decay
    double sum = 0.0;
    int n = 0;
    for (int ring = 1; ring <= R.max_used; ++ring) {
      if (R.total[ring] > 0) {
        sum += R.green[ring] / R.total[ring];
        ++n;
      }
    }
    return n > 0 ? sum / n : NA_REAL;
  }
  double num = 0.0, den = 0.0;
  for (int ring = 1; ring <= R.max_used; ++ring) {
    if (R.total[ring] > 0) {
      num += (R.green[ring] / R.total[ring]) * weight[ring];
      den += weight[ring];
    }
  }
  return num / den;
}

// Shared implementation of VGVI_cpp() and VGVI_rings_cpp().
struct VgviResult {
  std::vector<double> value;
  // only filled if rings are requested
  std::vector<std::vector<int>> ring, total;
  std::vector<std::vector<double>> green;
};

// Everything the parallel loop needs (plain C++ objects and raw pointers only).
struct VgviJob {
  const RasterInfo *dsm_ras;
  const Observers *obs;
  const cgei::LosTable *T;
  const cgei::VisMask *M;
  const cgei::BlockMax *bm;
  const GreenspaceMap *gmap;
  const double *dsm_v;           // DSM values (double; observer cells, bounds)
  const double *gs_v;
  const double *h0_v;
  const std::vector<int> *ring_of;   // ring of (|dr|, |dc|), (r + 1) x (r + 1)
  const std::vector<double> *weight;
  int r, fun, nthreads, batch;
  bool early_stop, want_rings;
};

// The observers in Morton order, in batches of job.batch: sweep a batch, then
// read every observer's mask in raster order into its ring histogram. V is the
// storage type of the DSM used by the sweep (float if exact, else double).
template <class V>
void vgvi_loop(const VgviJob &job, const V *dsm_sweep, cgei::ParallelProgress &progress,
               VgviResult &res) {
  const RasterInfo &ras = *job.dsm_ras;
  const Observers &obs = *job.obs;
  const cgei::LosTable &T = *job.T;
  const cgei::VisMask &M = *job.M;
  const GreenspaceMap &gmap = *job.gmap;
  const std::vector<int> &ring_of = *job.ring_of;
  const int r = job.r, nrow = ras.nrow, ncol = ras.ncol;

  // Per-thread scratch, allocated outside the parallel region.
  std::vector<cgei::BatchScratch> scratch;
  std::vector<Rings> rings;
  scratch.reserve(job.nthreads);
  rings.reserve(job.nthreads);
  int max_ring = 1;
  for (const int ring : ring_of) max_ring = std::max(max_ring, ring);
  for (int t = 0; t < job.nthreads; ++t) {
    scratch.emplace_back(T, M, job.batch);
    rings.emplace_back(max_ring);
  }

  const int nvalid = static_cast<int>(obs.order.size());
  const int nbatch = (nvalid + job.batch - 1) / job.batch;

#ifdef _OPENMP
#pragma omp parallel for num_threads(job.nthreads) schedule(dynamic, 1)
#endif
  for (int b = 0; b < nbatch; ++b) {
    if (progress.aborted()) continue;
    const int tid = cgei::omp_thread_num();
    cgei::BatchScratch &S = scratch[tid];
    Rings &R = rings[tid];
    const int i0 = b * job.batch, i1 = std::min(nvalid, i0 + job.batch);

    // The observer cell is always visible; lines of sight only if the eye
    // level is above the surface.
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
      cgei::drain_mask(M, S.mask(i - i0), [&](const int dr, const int dc) {
        const int ring = ring_of[static_cast<std::size_t>(std::abs(dr)) * (r + 1) + std::abs(dc)];
        // dr / dc are only needed if the greenspace grid differs from the DSM grid
        const long long cell = cell0 + static_cast<long long>(dr) * ncol + dc;
        R.add(ring, gmap.same_grid ? gmap.value(job.gs_v, 0, 0, cell)
                                   : gmap.value(job.gs_v, row0 + dr, col0 + dc, cell));
      });
      res.value[k] = vgvi_index(R, job.fun, *job.weight);
      if (job.want_rings) {
        for (int ring = 1; ring <= R.max_used; ++ring) {
          if (R.total[ring] > 0) {
            res.ring[k].push_back(ring);
            res.total[k].push_back(R.total[ring]);
            res.green[k].push_back(R.green[ring]);
          }
        }
      }
      R.reset();
    }
    progress.tick(static_cast<unsigned long>(i1 - i0));
  }
}

VgviResult run_vgvi(const NumericVector &dsm, const NumericVector &dsm_values, const NumericVector &greenspace,
                    const NumericVector &greenspace_values, const IntegerVector &x0,
                    const IntegerVector &y0, const NumericVector &h0, const int radius,
                    const int fun, const double m, const double b, const int ncores,
                    const bool display_progress, const bool early_stop, const bool want_rings) {
  const RasterInfo dsm_ras(dsm);
  const RasterInfo gs_ras(greenspace);
  if (dsm_values.size() != static_cast<R_xlen_t>(dsm_ras.nrow) * dsm_ras.ncol)
    Rcpp::stop("dsm_values does not match the dimensions of dsm.");
  if (greenspace_values.size() != static_cast<R_xlen_t>(gs_ras.nrow) * gs_ras.ncol)
    Rcpp::stop("greenspace_values does not match the dimensions of greenspace.");
  if (radius < 1) Rcpp::stop("radius must be at least 1.");
  if (fun < 1 || fun > 3) Rcpp::stop("fun must be 1 (logit), 2 (exponential) or 3 (none).");

  const Observers obs = prepare_observers(x0, y0, dsm_ras);
  const int n = static_cast<int>(obs.cell.size());
  if (h0.size() != n) Rcpp::stop("h0 must have the same length as x0.");

  // Line-of-sight geometry
  const int r = static_cast<int>(std::round(radius / dsm_ras.res));
  if (r > cgei::kMaxRadius) Rcpp::stop("max_distance is too large for the resolution of the DSM.");
  cgei::LosTable T(r);
  if (!T.bind(dsm_ras.ncol)) Rcpp::stop("Raster too large for the radius.");
  const cgei::VisMask M(r);

  // Distance ring of every cell of the neighbourhood: distance in map units,
  // rounded, at least 1 (by symmetry a function of |dr| and |dc|).
  std::vector<int> ring_of(static_cast<std::size_t>(r + 1) * (r + 1));
  int max_ring = 1;
  for (int a = 0; a <= r; ++a) {
    for (int c = 0; c <= r; ++c) {
      const double d = std::sqrt(static_cast<double>(a * a + c * c));
      const int ring = std::max(1, static_cast<int>(std::lround(dsm_ras.res * d)));
      ring_of[static_cast<std::size_t>(a) * (r + 1) + c] = ring;
      max_ring = std::max(max_ring, ring);
    }
  }

  // Decay weight of every ring: integral of the decay function over the ring's
  // normalised distance interval (computed once instead of once per observer).
  std::vector<double> weight(max_ring + 1, 0.0);
  if (fun != 3) {
    const double min_dist = 1 / static_cast<double>(radius);
    for (int ring = 1; ring <= max_ring; ++ring) {
      const double d = ring / static_cast<double>(radius);
      weight[ring] = integrate(d - min_dist, d, 300, fun, m, b);
    }
  }

  const GreenspaceMap gmap(dsm_ras, gs_ras);
  const double *dsm_v = dsm_values.begin();
  const int nthreads = cgei::resolve_threads(ncores);
  // Block maxima for the early-termination bounds (an empty grid if disabled)
  const cgei::BlockMax bm(dsm_v, early_stop ? dsm_ras.nrow : 0, early_stop ? dsm_ras.ncol : 0,
                          nthreads);

  VgviResult res;
  res.value.assign(n, NA_REAL);
  if (want_rings) {
    res.ring.resize(n);
    res.total.resize(n);
    res.green.resize(n);
  }

  VgviJob job;
  job.dsm_ras = &dsm_ras;
  job.obs = &obs;
  job.T = &T;
  job.M = &M;
  job.bm = &bm;
  job.gmap = &gmap;
  job.dsm_v = dsm_v;
  job.gs_v = greenspace_values.begin();
  job.h0_v = h0.begin();
  job.ring_of = &ring_of;
  job.weight = &weight;
  job.r = r;
  job.fun = fun;
  job.nthreads = nthreads;
  job.batch = cgei::batch_size(M, static_cast<int>(obs.order.size()), nthreads);
  job.early_stop = early_stop;
  job.want_rings = want_rings;

  cgei::ParallelProgress progress(obs.order.size(), display_progress);
  const std::size_t ncell = static_cast<std::size_t>(dsm_ras.nrow) * dsm_ras.ncol;
  if (cgei::float_exact(dsm_v, ncell, nthreads)) {
    std::vector<float> dsm_f(ncell);
    cgei::to_float(dsm_v, ncell, dsm_f.data(), nthreads);
    vgvi_loop<float>(job, dsm_f.data(), progress, res);
  } else {
    vgvi_loop<double>(job, dsm_v, progress, res);
  }
  progress.finish();  // throws if the user interrupted
  return res;
}

}  // namespace

// Viewshed Greenness Visibility Index for every observer.
// x0 / y0: 1-based column / row of the observers in dsm, h0: eye level,
// radius: maximum distance in map units, fun: 1 = logit, 2 = exponential,
// 3 = no decay. Observers outside the raster get NA.
// [[Rcpp::export]]
std::vector<double> VGVI_cpp(const Rcpp::NumericVector &dsm, const Rcpp::NumericVector &dsm_values,
                             const Rcpp::NumericVector &greenspace, const Rcpp::NumericVector &greenspace_values,
                             const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                             const Rcpp::NumericVector &h0, const int radius,
                             const int fun, const double m, const double b,
                             const int ncores = 1, const bool display_progress = false,
                             const bool early_stop = true) {
  return run_vgvi(dsm, dsm_values, greenspace, greenspace_values, x0, y0, h0, radius, fun,
                  m, b, ncores, display_progress, early_stop, false).value;
}

// Diagnostic variant of VGVI_cpp() used by the tests: the distance-ring
// histogram (non-empty rings only) behind every observer's VGVI value.
// [[Rcpp::export]]
Rcpp::List VGVI_rings_cpp(const Rcpp::NumericVector &dsm, const Rcpp::NumericVector &dsm_values,
                          const Rcpp::NumericVector &greenspace, const Rcpp::NumericVector &greenspace_values,
                          const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                          const Rcpp::NumericVector &h0, const int radius,
                          const int ncores = 1, const bool early_stop = true) {
  VgviResult res = run_vgvi(dsm, dsm_values, greenspace, greenspace_values, x0, y0, h0,
                            radius, 3, 1.0, 1.0, ncores, false, early_stop, true);
  const int n = static_cast<int>(res.value.size());
  Rcpp::List out(n);
  for (int k = 0; k < n; ++k) {
    out[k] = Rcpp::DataFrame::create(Rcpp::Named("ring") = Rcpp::wrap(res.ring[k]),
                                     Rcpp::Named("n_visible") = Rcpp::wrap(res.total[k]),
                                     Rcpp::Named("green") = Rcpp::wrap(res.green[k]));
  }
  return out;
}

// OpenMP availability of this build (used by the R wrappers and the tests).
// [[Rcpp::export]]
Rcpp::List cgei_openmp_info() {
#ifdef _OPENMP
  return Rcpp::List::create(Rcpp::Named("openmp") = true,
                            Rcpp::Named("max_threads") = omp_get_max_threads(),
                            Rcpp::Named("num_procs") = omp_get_num_procs());
#else
  return Rcpp::List::create(Rcpp::Named("openmp") = false,
                            Rcpp::Named("max_threads") = 1,
                            Rcpp::Named("num_procs") = 1);
#endif
}
