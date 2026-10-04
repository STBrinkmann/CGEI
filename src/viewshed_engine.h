// Viewshed engine shared by VGVI and VVI.
//
// Plain C++ (no Rcpp / R API): everything in here is safe to call from OpenMP
// worker threads, and the header compiles stand-alone for the sanitizer tests
// in dev/cpp-tests.
//
// Algorithm (unchanged from the original implementation): 8 * r Bresenham lines
// of sight radiate from the observer; walking a line outwards, a cell is
// visible if the tangent (h - h0) / d is strictly larger than the largest
// tangent seen so far on that line (the horizon). Consecutive lines share their
// first cells, so the horizon of the shared prefix is reused.
//
// Bugs fixed with respect to the original loop:
//  * the horizon of a shared prefix is reused for prefixes of length 1 as well
//    (the original used `k_i > 1` and restarted those lines at -9999, ignoring
//    an obstacle directly next to the observer);
//  * NA cells no longer skip the horizon book-keeping (the original left a
//    stale horizon of an unrelated line behind);
//  * cells outside the raster are detected from their row / column (the
//    original missed column wrap-arounds whenever the raster had <= 2r columns).
//
// Performance:
//  * all geometry (offsets, distances, prefix lengths) is precomputed once;
//  * no allocation per observer; first-visit de-duplication uses a reusable
//    byte mask that is reset through the list of touched cells;
//  * exact early termination: a line is abandoned as soon as no remaining cell
//    can beat the current horizon, using a cheap per-quadrant upper bound of
//    the DSM (block maxima). Results are identical with and without it.

#ifndef CGEI_VIEWSHED_ENGINE_H
#define CGEI_VIEWSHED_ENGINE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "los_geometry.h"

namespace cgei {

constexpr double kNoHorizon = -9999.0;  // initial horizon, as in the original code

// Precomputed line-of-sight table for a radius of r cells.
struct LosTable {
  int r = 0;
  int nref = 1;     // number of cells of the (2r+1) x (2r+1) reference grid
  int nlines = 0;

  // Lines in CSR layout: steps of line l are line_beg[l] .. line_beg[l + 1] - 1
  std::vector<int> line_beg;
  std::vector<int> start;           // first step not shared with the previous line
  std::vector<int> quadrant;        // quadrant (0..3) used for the early-termination bound
  std::vector<double> dist_last;    // distance of the last step of the line
  std::vector<unsigned char> monotone;  // distance strictly increasing along the line?

  // Per step
  std::vector<int> dr, dc;    // row / column offset from the observer
  std::vector<int> ref;       // reference-grid id (dr + r) * (2r + 1) + (dc + r)
  std::vector<double> dist;   // sqrt(dr^2 + dc^2) in cells (same expression as the original)

  // Hot data of the sweep, packed per step (set by bind())
  struct Step {
    double dist;   // as above
    int off;       // linear offset dr * ncol + dc
    int ref;       // as above
  };
  std::vector<Step> step;

  // Bounding box (in offsets) of all steps of each quadrant
  int qbox[4][4];             // [q][0..3] = min_dr, max_dr, min_dc, max_dc

  // All reference cells reached by any line (sorted by ref), used for VVI
  std::vector<int> cover_dr, cover_dc;

  double max_dist = 0.0;
  int max_len = 0;

  explicit LosTable(const int radius_cells) { build(radius_cells); }

  void build(const int radius_cells) {
    r = radius_cells > 0 ? radius_cells : 0;
    const int nc_ref = 2 * r + 1;
    nref = nc_ref * nc_ref;
    nlines = 8 * r;
    line_beg.assign(1, 0);
    start.clear(); quadrant.clear(); dist_last.clear(); monotone.clear();
    dr.clear(); dc.clear(); ref.clear(); dist.clear(); step.clear();
    cover_dr.clear(); cover_dc.clear();
    max_dist = 0.0;
    max_len = 0;
    for (int q = 0; q < 4; ++q) {
      qbox[q][0] = qbox[q][2] = std::numeric_limits<int>::max();
      qbox[q][1] = qbox[q][3] = std::numeric_limits<int>::min();
    }
    if (r == 0) return;

    const std::vector<int> los = los_reference(r, r, r, nc_ref);
    const std::vector<int> shared = shared_los(r, los);
    std::vector<unsigned char> covered(static_cast<std::size_t>(nref), 0);

    for (int l = 0; l < nlines; ++l) {
      const int q = l / (2 * r);
      double prev = -1.0;
      bool mono = true;
      int len = 0;
      for (int j = 0; j < r; ++j) {
        const int cell = los[static_cast<std::size_t>(l) * r + j];
        if (cell == kNaInt) break;
        const int rr = cell / nc_ref, cc = cell - (cell / nc_ref) * nc_ref;
        const int ddr = rr - r, ddc = cc - r;
        const double d = std::sqrt(static_cast<double>(ddr * ddr + ddc * ddc));
        dr.push_back(ddr);
        dc.push_back(ddc);
        ref.push_back(cell);
        dist.push_back(d);
        if (!(d > prev)) mono = false;
        prev = d;
        if (d > max_dist) max_dist = d;
        covered[cell] = 1;
        qbox[q][0] = std::min(qbox[q][0], ddr);
        qbox[q][1] = std::max(qbox[q][1], ddr);
        qbox[q][2] = std::min(qbox[q][2], ddc);
        qbox[q][3] = std::max(qbox[q][3], ddc);
        ++len;
      }
      line_beg.push_back(static_cast<int>(dr.size()));
      // A shared prefix can never be longer than the line itself.
      start.push_back(std::min(shared[l], len));
      quadrant.push_back(q);
      dist_last.push_back(len > 0 ? prev : 0.0);
      monotone.push_back(mono ? 1 : 0);
      if (len > max_len) max_len = len;
    }
    for (int cell = 0; cell < nref; ++cell) {
      if (covered[cell]) {
        cover_dr.push_back(cell / nc_ref - r);
        cover_dc.push_back(cell % nc_ref - r);
      }
    }
  }

  // Bind the table to a raster with ncol columns (linear offsets).
  // Returns false if the offsets do not fit into an int.
  bool bind(const int ncol) {
    if (static_cast<double>(r) * ncol + r >= 2147483647.0) return false;
    step.resize(dr.size());
    for (std::size_t s = 0; s < dr.size(); ++s) {
      step[s].dist = dist[s];
      step[s].off = dr[s] * ncol + dc[s];
      step[s].ref = ref[s];
    }
    return true;
  }

  int ref_center() const { return r * (2 * r + 1) + r; }
};

// Block maxima of the DSM, used for an upper bound of the heights in a box.
struct BlockMax {
  static constexpr int B = 16;
  int nrow = 0, ncol = 0, nbr = 0, nbc = 0;
  std::vector<double> bmax;  // -inf for blocks without any non-NA cell

  BlockMax(const double *dsm, const int nrow_, const int ncol_, const int nthreads = 1) {
    nrow = nrow_;
    ncol = ncol_;
    nbr = (nrow + B - 1) / B;
    nbc = (ncol + B - 1) / B;
    bmax.assign(static_cast<std::size_t>(nbr) * nbc, -std::numeric_limits<double>::infinity());
    // One block row per task: every block is written by exactly one thread.
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (int br = 0; br < nbr; ++br) {
      double *bm = &bmax[static_cast<std::size_t>(br) * nbc];
      const int r1 = std::min(nrow, (br + 1) * B);
      for (int row = br * B; row < r1; ++row) {
        const double *v = dsm + static_cast<std::size_t>(row) * ncol;
        for (int bc = 0; bc < nbc; ++bc) {
          double m = bm[bc];
          const int c1 = std::min(ncol, (bc + 1) * B);
          for (int col = bc * B; col < c1; ++col) {
            if (v[col] > m) m = v[col];  // NaN never compares greater
          }
          bm[bc] = m;
        }
      }
    }
  }

  // Upper bound of the DSM in rows [r0, r1] x cols [c0, c1] (clipped to the raster).
  double query(int r0, int r1, int c0, int c1) const {
    if (r0 < 0) r0 = 0;
    if (c0 < 0) c0 = 0;
    if (r1 > nrow - 1) r1 = nrow - 1;
    if (c1 > ncol - 1) c1 = ncol - 1;
    double m = -std::numeric_limits<double>::infinity();
    if (r0 > r1 || c0 > c1) return m;
    for (int br = r0 / B; br <= r1 / B; ++br) {
      const double *bm = &bmax[static_cast<std::size_t>(br) * nbc];
      for (int bc = c0 / B; bc <= c1 / B; ++bc) {
        if (bm[bc] > m) m = bm[bc];
      }
    }
    return m;
  }
};

// Description of a raster grid (row 0 at the top / ymax).
struct GridInfo {
  int nrow = 0, ncol = 0;
  double xmin = 0, xmax = 0, ymin = 0, ymax = 0, res = 0;
  long long ncell() const { return static_cast<long long>(nrow) * ncol; }
};

// Per-thread scratch space of the sweep.
// Bytes of padding around per-thread objects that are stored next to each
// other (e.g. in a std::vector, one element per thread): their members are
// written in the innermost loops, and without padding the objects of
// neighbouring threads would share cache lines ("false sharing"). Two cache
// lines, as Intel CPUs prefetch pairs of lines.
constexpr std::size_t kCachePad = 128;

struct SweepScratch {
  char pad_front_[kCachePad];
  std::vector<double> horizon;        // horizon (max tangent) per step of the current line
  std::vector<unsigned char> mask;    // first-visit flags on the reference grid
  std::vector<int> touched;           // reference ids set in `mask` (for resetting)
  char pad_back_[kCachePad];

  explicit SweepScratch(const LosTable &T)
      : horizon(static_cast<std::size_t>(std::max(1, T.max_len)), kNoHorizon),
        mask(static_cast<std::size_t>(T.nref), 0) {
    touched.reserve(1024);
  }

  // Mark a reference cell as visited; returns true on the first visit.
  inline bool first_visit(const int ref) {
    if (mask[ref]) return false;
    mask[ref] = 1;
    touched.push_back(ref);
    return true;
  }

  inline void reset() {
    for (int ref : touched) mask[ref] = 0;
    touched.clear();
  }
};

// Upper bounds of (DSM - h0) for the four quadrants around an observer;
// -inf if the quadrant holds no valid cell.
inline void quadrant_bounds(const LosTable &T, const BlockMax &bm, const int row0,
                            const int col0, const double h0, double out[4]) {
  for (int q = 0; q < 4; ++q) {
    if (T.qbox[q][0] > T.qbox[q][1]) {  // empty quadrant
      out[q] = -std::numeric_limits<double>::infinity();
      continue;
    }
    const double hmax = bm.query(row0 + T.qbox[q][0], row0 + T.qbox[q][1],
                                 col0 + T.qbox[q][2], col0 + T.qbox[q][3]);
    out[q] = hmax - h0;
  }
}

// Distance threshold beyond which no cell of the current line can be visible:
// with a = upper bound of (h - h0) and horizon t, any cell at distance d has
// tangent <= a / d. For a > 0 that is <= t once d >= a / t (the factor
// 1 + 1e-12 makes the test conservative w.r.t. rounding, so the result is
// exactly the one without early termination). For a <= 0 the largest
// possible tangent of the remaining cells is a / d_last.
inline double stop_distance(const double a, const double horizon, const double dist_last) {
  const double inf = std::numeric_limits<double>::infinity();
  if (a > 0) {
    return (horizon > 0) ? (a / horizon) * (1.0 + 1e-12) : inf;
  }
  if (a <= 0) {  // false for NaN
    return (horizon >= a / dist_last) ? 0.0 : inf;
  }
  return inf;
}

// Sweep all lines of sight of one observer.
//
//  dsm      DSM values (row-major, NaN = NA)
//  cell0    linear index of the observer cell, row0 / col0 its row / column
//  h0       observer eye level (must be above the DSM at the observer cell;
//           the caller checks this)
//  qbound   quadrant bounds from quadrant_bounds() (ignored if !early_stop)
//  visit    called as visit(s, cell) whenever a cell is visible on a line (s is
//           the step index into the table, cell the linear raster index); a cell
//           can be reported by several lines, use first_visit() to de-duplicate
template <bool CheckBounds, class Visit>
inline void sweep_lines(const LosTable &T, const double *dsm, const int nrow,
                        const int ncol, const long long cell0, const int row0,
                        const int col0, const double h0, const double *qbound,
                        const bool early_stop, SweepScratch &S, Visit &&visit) {
  double *hz = S.horizon.data();
  int valid_upto = -1;  // hz[0..valid_upto] hold the horizon of the previous line
  const double inf = std::numeric_limits<double>::infinity();

  for (int l = 0; l < T.nlines; ++l) {
    const int beg = T.line_beg[l];
    const int end = T.line_beg[l + 1];
    const int k = T.start[l];

    double horizon;
    if (k == 0) {
      horizon = kNoHorizon;
    } else {
      if (k - 1 > valid_upto) {
        // The previous line stopped early, i.e. none of its remaining cells
        // (which include our shared prefix) could raise its horizon.
        const double fill = valid_upto >= 0 ? hz[valid_upto] : kNoHorizon;
        for (int t = valid_upto + 1; t < k; ++t) hz[t] = fill;
      }
      horizon = hz[k - 1];
    }

    const bool can_stop = early_stop && T.monotone[l];
    const double a = can_stop ? qbound[T.quadrant[l]] : 0.0;
    double dstop = can_stop ? stop_distance(a, horizon, T.dist_last[l]) : inf;

    const LosTable::Step *step = T.step.data();
    int j = k;
    int s = beg + k;
    for (; s < end; ++s, ++j) {
      const double d = step[s].dist;
      if (d >= dstop) break;
      if (CheckBounds) {
        const int rr = row0 + T.dr[s];
        const int cc = col0 + T.dc[s];
        if (static_cast<unsigned>(rr) >= static_cast<unsigned>(nrow) ||
            static_cast<unsigned>(cc) >= static_cast<unsigned>(ncol)) {
          hz[j] = horizon;
          continue;
        }
      }
      const long long cell = cell0 + step[s].off;
      const double h = dsm[cell];
      if (!std::isnan(h)) {
        const double tangent = (h - h0) / d;
        if (tangent > horizon) {
          horizon = tangent;
          visit(s, cell);
          if (can_stop) dstop = stop_distance(a, horizon, T.dist_last[l]);
        }
      }
      hz[j] = horizon;
    }
    valid_upto = (s < end) ? j - 1 : (end - beg) - 1;
  }
}

// Potential viewshed (VVI): the cells reached by any line of sight (static
// geometry, LosTable::cover_*) as runs of consecutive columns per row offset,
// in row-major order. The observer cell itself is not part of the runs.
struct CoverRuns {
  struct Run {
    int dr, lo, hi;
  };
  std::vector<Run> runs;

  explicit CoverRuns(const LosTable &T) {
    for (std::size_t c = 0; c < T.cover_dr.size(); ++c) {
      const int dr = T.cover_dr[c], dc = T.cover_dc[c];
      if (dr == 0 && dc == 0) continue;
      if (!runs.empty() && runs.back().dr == dr && runs.back().hi + 1 == dc) {
        runs.back().hi = dc;
      } else {
        runs.push_back({dr, dc, dc});
      }
    }
  }

  // Calls f(row, c0, c1) for the part of every run that lies inside an
  // nrow x ncol raster, for an observer at (row0, col0); row-major order.
  template <class F>
  void for_each(const int row0, const int col0, const int nrow, const int ncol, F f) const {
    for (const Run &run : runs) {
      const int rr = row0 + run.dr;
      if (rr < 0 || rr >= nrow) continue;
      const int c0 = std::max(col0 + run.lo, 0);
      const int c1 = std::min(col0 + run.hi, ncol - 1);
      if (c0 <= c1) f(rr, c0, c1);
    }
  }
};

// Number of valid (non-NaN) cells in a part of a row, from row-wise prefix
// counts (only stored if the raster has NaN cells).
struct ValidPrefix {
  int ncol = 0;
  bool all_valid = true;
  std::vector<int> prefix;  // nrow x (ncol + 1): valid cells left of column c

  ValidPrefix(const double *v, const int nrow, const int ncol_, const int nthreads) : ncol(ncol_) {
    const long long n = static_cast<long long>(nrow) * ncol;
    int has_nan = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) reduction(| : has_nan) schedule(static)
#endif
    for (long long i = 0; i < n; ++i) has_nan |= std::isnan(v[i]) ? 1 : 0;
    (void)nthreads;
    all_valid = has_nan == 0;
    if (all_valid) return;
    prefix.assign(static_cast<std::size_t>(nrow) * (ncol + 1), 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (int row = 0; row < nrow; ++row) {
      const double *src = v + static_cast<std::size_t>(row) * ncol;
      int *dst = prefix.data() + static_cast<std::size_t>(row) * (ncol + 1);
      int s = 0;
      dst[0] = 0;
      for (int c = 0; c < ncol; ++c) {
        s += std::isnan(src[c]) ? 0 : 1;
        dst[c + 1] = s;
      }
    }
  }

  int count(const int row, const int c0, const int c1) const {
    if (all_valid) return c1 - c0 + 1;
    const int *p = prefix.data() + static_cast<std::size_t>(row) * (ncol + 1);
    return p[c1 + 1] - p[c0];
  }
};

// Size of the potential viewshed of an observer at (row0, col0): the covered
// cells inside the raster with a valid height, plus the observer cell.
inline int count_seen(const CoverRuns &cov, const ValidPrefix &valid, const int row0,
                      const int col0, const int nrow, const int ncol) {
  int n = 1;
  cov.for_each(row0, col0, nrow, ncol,
               [&](const int rr, const int c0, const int c1) { n += valid.count(rr, c0, c1); });
  return n;
}

// Per cell, the number of observers (rows[k], cols[k]) whose potential
// viewshed contains the cell (out: nrow * ncol values, overwritten), i.e. the
// counts of count_seen()'s cells over all observers. Uses difference arrays
// along the rows: O(observers x runs + cells).
inline void accumulate_seen(const CoverRuns &cov, const double *dsm, const int nrow, const int ncol,
                            const std::vector<int> &rows, const std::vector<int> &cols,
                            const int nthreads, int *out) {
  const std::size_t stride = static_cast<std::size_t>(ncol) + 1;
  std::vector<int> diff(static_cast<std::size_t>(nrow) * stride, 0);
  for (std::size_t k = 0; k < rows.size(); ++k) {
    cov.for_each(rows[k], cols[k], nrow, ncol, [&](const int rr, const int c0, const int c1) {
      int *d = diff.data() + static_cast<std::size_t>(rr) * stride;
      ++d[c0];
      --d[c1 + 1];
    });
  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int row = 0; row < nrow; ++row) {
    const int *d = diff.data() + static_cast<std::size_t>(row) * stride;
    const double *h = dsm + static_cast<std::size_t>(row) * ncol;
    int *o = out + static_cast<std::size_t>(row) * ncol;
    int s = 0;
    for (int c = 0; c < ncol; ++c) {
      s += d[c];
      o[c] = std::isnan(h[c]) ? 0 : s;
    }
  }
  (void)nthreads;
  for (std::size_t k = 0; k < rows.size(); ++k) {
    ++out[static_cast<std::size_t>(rows[k]) * ncol + cols[k]];
  }
}

// Morton (Z-order) key for cache-friendly processing order of observers.
inline std::uint64_t morton_key(const std::uint32_t row, const std::uint32_t col) {
  auto spread = [](std::uint64_t v) {
    v &= 0xFFFFFFFFULL;
    v = (v | (v << 16)) & 0x0000FFFF0000FFFFULL;
    v = (v | (v << 8)) & 0x00FF00FF00FF00FFULL;
    v = (v | (v << 4)) & 0x0F0F0F0F0F0F0F0FULL;
    v = (v | (v << 2)) & 0x3333333333333333ULL;
    v = (v | (v << 1)) & 0x5555555555555555ULL;
    return v;
  };
  return spread(col) | (spread(row) << 1);
}

}  // namespace cgei

#endif
