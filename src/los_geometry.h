// Line-of-sight (LoS) reference geometry.
//
// Plain C++ (no Rcpp / R API) so that it can be used from OpenMP worker threads
// and compiled stand-alone for sanitizer tests (dev/cpp-tests).
//
// The geometry is a faithful port of the original bresenham_map() /
// LoS_reference() / shared_LoS() code and produces bit-identical output.

#ifndef CGEI_LOS_GEOMETRY_H
#define CGEI_LOS_GEOMETRY_H

#include <climits>
#include <cstddef>
#include <vector>

namespace cgei {

// R's NA_integer_
constexpr int kNaInt = INT_MIN;

// Bresenham lines of the first octant: from (x0, y0) towards (x0 + i, y0 + radius)
// for i = 0..radius, truncated at the circle of the given radius (in cells).
// Returns a row-major (radius + 1) x radius matrix of cell ids y * nc + x,
// padded with kNaInt.
inline std::vector<int> bresenham_lines(const int x0, const int y0,
                                        const int radius, const int nc) {
  std::vector<int> out(static_cast<std::size_t>(radius + 1) * radius, kNaInt);
  const int Dy = radius;
  const int Sy = (Dy > 0) - (Dy < 0);
  for (int i = 0; i <= radius; ++i) {
    const int Dx = i;
    const int Sx = (Dx > 0) - (Dx < 0);
    int R = radius / 2;  // initial remainder (integer division, as in the original)
    int X = x0, Y = y0;
    Y += Sy;
    R += Dx;
    if (R >= Dy) {
      X += Sx;
      R -= Dy;
    }
    int c = 0;
    while (c < radius && ((x0 - X) * (x0 - X) + (y0 - Y) * (y0 - Y)) <= radius * radius) {
      out[static_cast<std::size_t>(i) * radius + c] = Y * nc + X;
      Y += Sy;
      R += Dx;
      if (R >= Dy) {
        X += Sx;
        R -= Dy;
      }
      ++c;
    }
  }
  return out;
}

// Reference LoS paths for a (2r+1) x (2r+1) reference grid centred at
// (x0_ref, y0_ref): 8 * r lines with r steps each, row-major [line * r + step],
// holding 0-based reference cell ids (row * nc_ref + col), padded with kNaInt.
// Lines are ordered clockwise; lines 0..r form the first octant and the other
// seven octants are rotations / reflections of it.
inline std::vector<int> los_reference(const int x0_ref, const int y0_ref,
                                      const int r, const int nc_ref) {
  std::vector<int> los(static_cast<std::size_t>(8) * r * r, 0);
  if (r <= 0) return los;
  const int l = r + 1;
  const std::vector<int> bh = bresenham_lines(x0_ref, y0_ref, r, nc_ref);
  auto set = [&](const int line, const int step, const int value) {
    los[static_cast<std::size_t>(line) * r + step] = value;
  };
  for (int i = 0; i < l; ++i) {
    for (int j = 0; j < r; ++j) {
      const int cell = bh[static_cast<std::size_t>(i) * r + j];
      if (cell == kNaInt) {
        set(0 * r + i, j, kNaInt);
        set(2 * r + i, j, kNaInt);
        set(4 * r + i, j, kNaInt);
        set(6 * r + i, j, kNaInt);
        if (i != 0 && i != (l - 1)) {
          set(2 * r - i, j, kNaInt);
          set(4 * r - i, j, kNaInt);
          set(6 * r - i, j, kNaInt);
          set(8 * r - i, j, kNaInt);
        }
      } else {
        const int row = cell / nc_ref;
        const int col = cell - row * nc_ref;
        const int x = col - x0_ref;
        const int y = row - y0_ref;
        set(0 * r + i, j, (y + y0_ref) * nc_ref + (x + x0_ref));
        set(2 * r + i, j, (-x + y0_ref) * nc_ref + (y + x0_ref));
        set(4 * r + i, j, (-y + y0_ref) * nc_ref + (-x + x0_ref));
        set(6 * r + i, j, (x + y0_ref) * nc_ref + (-y + x0_ref));
        if (i != 0 && i != (l - 1)) {
          set(2 * r - i, j, (x + y0_ref) * nc_ref + (y + x0_ref));
          set(4 * r - i, j, (-y + y0_ref) * nc_ref + (x + x0_ref));
          set(6 * r - i, j, (-x + y0_ref) * nc_ref + (-y + x0_ref));
          set(8 * r - i, j, (y + y0_ref) * nc_ref + (-x + x0_ref));
        }
      }
    }
  }
  return los;
}

// For every line, the first step at which it differs from the previous line
// (0 for the first line, and for a line identical to its predecessor).
// Steps before that index are shared with the previous line, so their horizon
// (maximum tangent) can be reused.
inline std::vector<int> shared_los(const int r, const std::vector<int> &los) {
  std::vector<int> out(static_cast<std::size_t>(8) * (r > 0 ? r : 0), 0);
  for (int i = 1; i < 8 * r; ++i) {
    for (int j = 0; j < r; ++j) {
      if (los[static_cast<std::size_t>(i) * r + j] !=
          los[static_cast<std::size_t>(i - 1) * r + j]) {
        out[i] = j;
        break;
      }
    }
  }
  return out;
}

}  // namespace cgei

#endif
