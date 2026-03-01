# Weighted excess kurtosis

Computes the population weighted excess kurtosis:
`sum(w_hat * ((x - mu) / sigma)^4) - 3`.

## Usage

``` r
weighted_kurtosis(x, w, na.rm = FALSE)
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
