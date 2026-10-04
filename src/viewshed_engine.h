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
//  * all geometry (offsets, distances, prefix lengths) is precomputed once and
//    packed into 16 bytes per step;
//  * exact early termination: a line is abandoned as soon as no remaining cell
//    can beat the current horizon, using a cheap per-quadrant upper bound of
//    the DSM (block maxima). Results are identical with and without it;
//  * visible cells are only flagged in a per-observer bit mask of the
//    (2r+1) x (2r+1) neighbourhood (one OR per visible step, which also
//    de-duplicates cells seen by several lines). The callers read the mask
//    afterwards in raster order (drain_mask()), i.e. with sequential access to
//    their own rasters instead of a random access per visible step;
//  * observers are swept in batches of nearby observers (Morton order), line by
//    line: the steps of a line are loaded once per batch and the DSM cells of
//    neighbouring observers are mostly the same, so both stay in the caches;
//  * the DSM is read as float if all its values are exactly representable as
//    float (e.g. read from a Float32 GeoTIFF): half the memory traffic, same
//    values and therefore identical results.
//
// The order in which the lines of an observer are walked, and the order of the
// observers, do not change the result: the visibility of a cell on a line only
// depends on the cells before it on that line, and the mask is the union over
// all lines.

#ifndef CGEI_VIEWSHED_ENGINE_H
#define CGEI_VIEWSHED_ENGINE_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "los_geometry.h"

namespace cgei {

constexpr double kNoHorizon = -9999.0;  // initial horizon, as in the original code

// Largest supported radius in cells (row / column offsets are stored as int16).
constexpr int kMaxRadius = 32767;

// Precomputed line-of-sight table for a radius of r cells.
struct LosTable {
  int r = 0;
  int nlines = 0;
  bool too_large = false;  // radius above kMaxRadius (empty table)

  // Lines in CSR layout: steps of line l are line_beg[l] .. line_beg[l + 1] - 1
  std::vector<int> line_beg;
  std::vector<int> start;           // first step not shared with the previous line
  std::vector<int> quadrant;        // quadrant (0..3) used for the early-termination bound
  std::vector<double> dist_last;    // distance of the last step of the line
  std::vector<unsigned char> monotone;  // distance strictly increasing along the line?

  // Hot data of the sweep, one record per step
  struct Step {
    double dist;          // sqrt(dr^2 + dc^2) in cells (same expression as the original)
    int off;              // linear offset dr * ncol + dc (set by bind())
    std::int16_t dr, dc;  // row / column offset from the observer
  };
  std::vector<Step> step;

  // Bounding box (in offsets) of all steps of each quadrant
  int qbox[4][4];             // [q][0..3] = min_dr, max_dr, min_dc, max_dc

  // All reference cells reached by any line (row-major), used for VVI
  std::vector<int> cover_dr, cover_dc;

  int max_len = 0;

  explicit LosTable(const int radius_cells) { build(radius_cells); }

  void build(const int radius_cells) {
    r = radius_cells > 0 ? radius_cells : 0;
    too_large = r > kMaxRadius;
    if (too_large) r = 0;
    const int nc_ref = 2 * r + 1;
    nlines = 8 * r;
    line_beg.assign(1, 0);
    start.clear(); quadrant.clear(); dist_last.clear(); monotone.clear();
    step.clear(); cover_dr.clear(); cover_dc.clear();
    max_len = 0;
    for (int q = 0; q < 4; ++q) {
      qbox[q][0] = qbox[q][2] = std::numeric_limits<int>::max();
      qbox[q][1] = qbox[q][3] = std::numeric_limits<int>::min();
    }
    if (r == 0) return;

    const std::vector<int> los = los_reference(r, r, r, nc_ref);
    const std::vector<int> shared = shared_los(r, los);
    std::vector<unsigned char> covered(static_cast<std::size_t>(nc_ref) * nc_ref, 0);

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
        step.push_back(Step{d, 0, static_cast<std::int16_t>(ddr), static_cast<std::int16_t>(ddc)});
        if (!(d > prev)) mono = false;
        prev = d;
        covered[cell] = 1;
        qbox[q][0] = std::min(qbox[q][0], ddr);
        qbox[q][1] = std::max(qbox[q][1], ddr);
        qbox[q][2] = std::min(qbox[q][2], ddc);
        qbox[q][3] = std::max(qbox[q][3], ddc);
        ++len;
      }
      line_beg.push_back(static_cast<int>(step.size()));
      // A shared prefix can never be longer than the line itself.
      start.push_back(std::min(shared[l], len));
      quadrant.push_back(q);
      dist_last.push_back(len > 0 ? prev : 0.0);
      monotone.push_back(mono ? 1 : 0);
      if (len > max_len) max_len = len;
    }
    for (int cell = 0; cell < nc_ref * nc_ref; ++cell) {
      if (covered[cell]) {
        cover_dr.push_back(cell / nc_ref - r);
        cover_dc.push_back(cell % nc_ref - r);
      }
    }
  }

  // Bind the table to a raster with ncol columns (linear offsets).
  // Returns false if the radius is too large or the offsets do not fit into an int.
  bool bind(const int ncol) {
    if (too_large) return false;
    if (static_cast<double>(r) * ncol + r >= 2147483647.0) return false;
    for (Step &s : step) s.off = s.dr * ncol + s.dc;
    return true;
  }
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
    (void)nthreads;
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

// Bytes of padding around per-thread objects that are stored next to each
// other (e.g. in a std::vector, one element per thread): their members are
// written in the innermost loops, and without padding the objects of
// neighbouring threads would share cache lines ("false sharing"). Two cache
// lines, as Intel CPUs prefetch pairs of lines.
constexpr std::size_t kCachePad = 128;

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

// ---------------------------------------------------------------------------
// DSM storage
// ---------------------------------------------------------------------------

// True if every value is NA, infinite or exactly representable as float, i.e.
// if a float copy holds exactly the same values.
inline bool float_exact(const double *v, const std::size_t n, const int nthreads) {
  int bad = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) reduction(| : bad) schedule(static)
#endif
  for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n); ++i) {
    const double x = v[i];
    if (std::isnan(x) || std::isinf(x)) continue;
    // out-of-range conversions to float are undefined: test the range first
    bad |= (std::fabs(x) > static_cast<double>(FLT_MAX) ||
            static_cast<double>(static_cast<float>(x)) != x) ? 1 : 0;
  }
  (void)nthreads;
  return bad == 0;
}

inline void to_float(const double *v, const std::size_t n, float *out, const int nthreads) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n); ++i) out[i] = static_cast<float>(v[i]);
  (void)nthreads;
}

// ---------------------------------------------------------------------------
// Visibility masks
// ---------------------------------------------------------------------------

// Bit mask of the (2r+1) x (2r+1) neighbourhood of an observer; each row is
// padded to whole 64-bit words, so that the mask can be read row by row.
struct VisMask {
  int r = 0;
  int W = 0;            // 64-bit words per row
  int rowlen = 0;       // bits per row (W * 64)
  int center = 0;       // bit of the observer cell
  std::size_t words = 0;

  explicit VisMask(const int r_) : r(r_), W((2 * r_ + 1 + 63) / 64) {
    rowlen = W * 64;
    center = r * rowlen + r;
    words = static_cast<std::size_t>(2 * r + 1) * W;
  }
  inline int bit(const int dr, const int dc) const { return center + dr * rowlen + dc; }
};

inline void mask_set(std::uint64_t *m, const int bit) {
  m[bit >> 6] |= std::uint64_t(1) << (bit & 63);
}

inline int popcount64(std::uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_popcountll(x);
#else
  int n = 0;
  for (; x; x &= x - 1) ++n;
  return n;
#endif
}

inline int ctz64(const std::uint64_t x) {  // x != 0
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctzll(x);
#else
  int n = 0;
  while (!((x >> n) & 1)) ++n;
  return n;
#endif
}

// Calls f(dr, dc) for every flagged cell in raster (row-major) order and
// clears the mask.
template <class F>
inline void drain_mask(const VisMask &M, std::uint64_t *m, F &&f) {
  for (int row = 0; row < 2 * M.r + 1; ++row) {
    std::uint64_t *w = m + static_cast<std::size_t>(row) * M.W;
    for (int i = 0; i < M.W; ++i) {
      std::uint64_t x = w[i];
      if (!x) continue;
      w[i] = 0;
      do {
        const int b = ctz64(x);
        x &= x - 1;
        f(row - M.r, i * 64 + b - M.r);
      } while (x);
    }
  }
}

// Number of flagged cells; clears the mask.
inline int count_mask(const VisMask &M, std::uint64_t *m) {
  int n = 0;
  for (std::size_t i = 0; i < M.words; ++i) {
    if (m[i]) {
      n += popcount64(m[i]);
      m[i] = 0;
    }
  }
  return n;
}

// ---------------------------------------------------------------------------
// Batched sweep
// ---------------------------------------------------------------------------

// One observer of a batch.
struct Viewer {
  long long cell0;        // linear index of the observer cell
  int row0, col0;
  double h0;              // eye level (above the DSM at the observer cell)
  double qbound[4];       // quadrant_bounds() (unused without early termination)
  bool interior;          // the whole circle lies inside the raster
  std::uint64_t *mask;    // VisMask::words words; visible cells are OR-ed in
  double *hz;             // LosTable::max_len values: horizon per step
  int valid_upto;         // hz[0..valid_upto] hold the horizon of the previous line
};

// Walk one line of sight of one observer, starting at step k (the steps
// before are shared with the previous line).
template <bool CheckBounds, class V>
inline void sweep_line(const LosTable::Step *step, const int beg, const int end, const int k,
                       const bool can_stop, const double a, const double dist_last,
                       const VisMask &M, const V *dsm, const int nrow, const int ncol,
                       Viewer &o) {
  double *hz = o.hz;
  double horizon;
  if (k == 0) {
    horizon = kNoHorizon;
  } else {
    if (k - 1 > o.valid_upto) {
      // The previous line stopped early, i.e. none of its remaining cells
      // (which include our shared prefix) could raise its horizon.
      const double fill = o.valid_upto >= 0 ? hz[o.valid_upto] : kNoHorizon;
      for (int t = o.valid_upto + 1; t < k; ++t) hz[t] = fill;
    }
    horizon = hz[k - 1];
  }
  double dstop = can_stop ? stop_distance(a, horizon, dist_last)
                          : std::numeric_limits<double>::infinity();
  const V *base = dsm + o.cell0;
  const double h0 = o.h0;
  std::uint64_t *mask = o.mask;
  int j = k;
  int s = beg + k;
  for (; s < end; ++s, ++j) {
    const double d = step[s].dist;
    if (d >= dstop) break;
    if (CheckBounds) {
      const int rr = o.row0 + step[s].dr;
      const int cc = o.col0 + step[s].dc;
      if (static_cast<unsigned>(rr) >= static_cast<unsigned>(nrow) ||
          static_cast<unsigned>(cc) >= static_cast<unsigned>(ncol)) {
        hz[j] = horizon;
        continue;
      }
    }
    const double h = static_cast<double>(base[step[s].off]);
    if (!std::isnan(h)) {
      const double tangent = (h - h0) / d;
      if (tangent > horizon) {
        horizon = tangent;
        mask_set(mask, M.bit(step[s].dr, step[s].dc));
        if (can_stop) dstop = stop_distance(a, horizon, dist_last);
      }
    }
    hz[j] = horizon;
  }
  o.valid_upto = (s < end) ? j - 1 : (end - beg) - 1;
}

// Sweep all lines of sight of n observers, line by line (the observers should
// be close to each other, e.g. consecutive in Morton order). Visible cells are
// OR-ed into the observers' masks; the observer cells themselves are not set.
//
//  dsm   DSM values (row-major, NaN = NA), float or double
template <class V>
inline void sweep_batch(const LosTable &T, const VisMask &M, const V *dsm, const int nrow,
                        const int ncol, const bool early_stop, Viewer *vw, const int n) {
  for (int a = 0; a < n; ++a) vw[a].valid_upto = -1;
  const LosTable::Step *step = T.step.data();
  for (int l = 0; l < T.nlines; ++l) {
    const int beg = T.line_beg[l];
    const int end = T.line_beg[l + 1];
    const int k = T.start[l];
    const bool can_stop = early_stop && T.monotone[l];
    const int q = T.quadrant[l];
    const double dl = T.dist_last[l];
    for (int a = 0; a < n; ++a) {
      Viewer &o = vw[a];
      const double bound = can_stop ? o.qbound[q] : 0.0;
      if (o.interior) {
        sweep_line<false>(step, beg, end, k, can_stop, bound, dl, M, dsm, nrow, ncol, o);
      } else {
        sweep_line<true>(step, beg, end, k, can_stop, bound, dl, M, dsm, nrow, ncol, o);
      }
    }
  }
}

// Number of observers per batch: at most 16, and at most about 8 MB of masks
// per thread for large radii.
inline int batch_size(const VisMask &M) {
  const std::size_t bytes = M.words * sizeof(std::uint64_t);
  const std::size_t budget = static_cast<std::size_t>(8) << 20;
  const std::size_t b = bytes > 0 ? budget / bytes : 16;
  return static_cast<int>(std::max<std::size_t>(1, std::min<std::size_t>(16, b)));
}

// Per-thread scratch of a batch: masks (zero between observers) and horizon
// arrays, padded against false sharing (see kCachePad).
struct BatchScratch {
  char pad_front_[kCachePad];
  int capacity;
  std::size_t words, hz_stride;
  std::vector<std::uint64_t> masks;
  std::vector<double> hz;
  std::vector<Viewer> viewers;
  char pad_back_[kCachePad];

  BatchScratch(const LosTable &T, const VisMask &M, const int capacity_)
      : capacity(capacity_), words(M.words),
        hz_stride(static_cast<std::size_t>(std::max(1, T.max_len)) + 8),
        masks(static_cast<std::size_t>(capacity_) * M.words, 0),
        hz(static_cast<std::size_t>(capacity_) * hz_stride, kNoHorizon) {
    viewers.reserve(capacity_);
  }
  std::uint64_t *mask(const int i) { return masks.data() + static_cast<std::size_t>(i) * words; }
  double *horizon(const int i) { return hz.data() + static_cast<std::size_t>(i) * hz_stride; }
};

// ---------------------------------------------------------------------------
// Potential viewshed (VVI)
// ---------------------------------------------------------------------------

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
