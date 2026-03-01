# Weighted quantiles

Computes weighted quantiles using a type-7 analog: observations are
placed at the midpoints of their weight intervals in the cumulative
weight distribution, rescaled to \[0, 1\], and quantiles are obtained by
linear interpolation. For equal weights this matches
`quantile(x, type = 7)`.

## Usage

``` r
weighted_quantile(x, w, probs = seq(0, 1, 0.25), na.rm = FALSE)
```

## Arguments

- x:

  A numeric vector of observations.

- w:

  A numeric vector of non-negative weights (need not sum to 1).

- probs:

  A numeric vector of probabilities in \[0, 1\].

- na.rm:

  Logical. If `TRUE`, paired `NA`s in `x` and `w` are removed.

## Value

A numeric vector of the same length as `probs`.
