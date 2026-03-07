# wstats

<!-- badges: start -->
[![R-CMD-check](https://github.com/jgaeb/wstats/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jgaeb/wstats/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/jgaeb/wstats/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/jgaeb/wstats/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

Base R provides `weighted.mean()` but nothing else. **wstats** fills the gap
with weighted versions of the other common descriptive statistics.

## Installation

```r
# install.packages("devtools")
devtools::install_github("jgaeb/wstats")
```

## Usage

```r
library(wstats)

x <- c(1, 2, 3, 4, 5)
w <- c(0.5, 1.0, 2.0, 1.0, 0.5)   # unnormalised importance weights

weighted_var(x, w)
weighted_sd(x, w)
weighted_quantile(x, w, probs = c(0.25, 0.5, 0.75))
weighted_median(x, w)
weighted_mad(x, w)
weighted_skewness(x, w)
weighted_kurtosis(x, w)

y <- c(2, 3, 1, 5, 4)
weighted_cov(x, y, w)
weighted_cor(x, y, w)
```

## Convention

All functions currently use a **population formula** — weights are treated as
probability masses of a discrete distribution, not survey sampling weights.
Concretely, variance is `Σ(ŵᵢ (xᵢ − μ)²)` where `ŵᵢ = wᵢ / Σwⱼ`, with no
Bessel correction. This is appropriate for importance sampling and the Bayesian
bootstrap; support for survey/frequency weights may be added in a future
release.

Computationally intensive routines are implemented in C++ via
[cpp11](https://cpp11.r-lib.org/).
