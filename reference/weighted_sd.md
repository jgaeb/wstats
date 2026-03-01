# Weighted standard deviation

Square root of
[`weighted_var()`](https://jgaeb.github.io/wstats/reference/weighted_var.md).

## Usage

``` r
weighted_sd(x, w, na.rm = FALSE)
```

## Arguments

- x:

  A numeric vector of observations.

- w:

  A numeric vector of non-negative weights (need not sum to 1).

- na.rm:

  Logical. If `TRUE`, paired `NA`s in `x` and `w` are removed before
  computation. If `FALSE` (default) and any `NA` is present, `NA` is
  returned.

## Value

A single numeric value.
