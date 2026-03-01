# Weighted correlation

Weighted correlation coefficient, computed as
`weighted_cov(x, y, w) / (weighted_sd(x, w) * weighted_sd(y, w))`.

## Usage

``` r
weighted_cor(x, y, w, na.rm = FALSE)
```

## Arguments

- x:

  A numeric vector of observations.

- y:

  A numeric vector of observations (same length as `x`).

- w:

  A numeric vector of non-negative weights (need not sum to 1).

- na.rm:

  Logical. If `TRUE`, observations with any `NA` in `x`, `y`, or `w` are
  removed before computation.

## Value

A single numeric value in \[-1, 1\].
