# CGEI 0.3.0

## New Features

-   `gavi()` function has been added. This function calculates the Greenspace Availability IndexIndex (GAVI) after a `lacunarity()` analysis.

## Docs

-   The C++ `lacunarity()` function now uses the R4 raster class.
-   The parameters `cores` and `progress` were adjusted in a uniform style in the existing functions.
-   Tests for the `gavi()` function have been added.


# CGEI 0.2.1

## Bug fixes

-   `vvi(mode = c("VVI", "cumulative"), by_row = TRUE)` now works as expected. ([#17](https://github.com/STBrinkmann/CGEI/issues/17))
-   `vvi()` manual has been adjusted. ([#23](https://github.com/STBrinkmann/CGEI/issues/23))

# CGEI 0.2.0

## New Features

-   `sf_interpolate_IDW()` function has been added. This is a highly efficient IDW method for converting an `sf` object to a `SpatRast`.

# CGEI 0.1.0

This is the first release of CGEI.

## New Features

-   `vvi()` function now accepts `mode = "cumulative"` and `by_row` arguments. ([#17](https://github.com/STBrinkmann/CGEI/issues/17), @hansvancalster)
