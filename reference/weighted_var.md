# Weighted variance

Computes the population (importance-weight) weighted variance
`sum(w_hat * (x - mu)^2)` where `w_hat = w / sum(w)` and
`mu = weighted.mean(x, w)`.

## Usage

``` r
weighted_var(x, w, na.rm = FALSE)
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
