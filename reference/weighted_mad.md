# Weighted median absolute deviation

Computes the weighted MAD as the weighted median of
`|x - median(x, w)|`, scaled by `constant` (default 1.4826 for
consistency with [`stats::mad()`](https://rdrr.io/r/stats/mad.html)
under normality).

## Usage

``` r
weighted_mad(x, w, na.rm = FALSE, constant = 1.4826)
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

- constant:

  Scale factor (default `1.4826`).

## Value

A single numeric value.
