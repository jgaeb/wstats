# ---- Input validation helpers ------------------------------------------------

.check_wx <- function(x, w, call = rlang::caller_env()) {
  if (!is.numeric(x))
    rlang::abort("`x` must be a numeric vector.", call = call)
  if (!is.numeric(w))
    rlang::abort("`w` must be a numeric vector.", call = call)
  if (length(x) != length(w))
    rlang::abort(
      glue::glue(
        "`x` and `w` must have the same length ({length(x)} vs {length(w)})."
      ),
      call = call
    )
  if (any(w < 0, na.rm = TRUE))
    rlang::abort("`w` must be non-negative.", call = call)
}

.check_wxy <- function(x, y, w, call = rlang::caller_env()) {
  if (!is.numeric(x))
    rlang::abort("`x` must be a numeric vector.", call = call)
  if (!is.numeric(y))
    rlang::abort("`y` must be a numeric vector.", call = call)
  if (!is.numeric(w))
    rlang::abort("`w` must be a numeric vector.", call = call)
  if (length(x) != length(w))
    rlang::abort(
      glue::glue(
        "`x` and `w` must have the same length ({length(x)} vs {length(w)})."
      ),
      call = call
    )
  if (length(y) != length(w))
    rlang::abort(
      glue::glue(
        "`y` and `w` must have the same length ({length(y)} vs {length(w)})."
      ),
      call = call
    )
  if (any(w < 0, na.rm = TRUE))
    rlang::abort("`w` must be non-negative.", call = call)
}

# NA handling: returns NULL if na.rm = FALSE and NAs are present (signal to
# return NA_real_); otherwise drops paired NAs from x and w and returns a list.
.drop_na_wx <- function(x, w, na.rm, call = rlang::caller_env()) {
  has_na <- anyNA(x) || anyNA(w)
  if (has_na) {
    if (!na.rm) return(NULL)
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }
  if (length(x) == 0L)
    rlang::abort("No non-NA observations after removing NAs.", call = call)
  list(x = x, w = w)
}

.drop_na_wxy <- function(x, y, w, na.rm, call = rlang::caller_env()) {
  has_na <- anyNA(x) || anyNA(y) || anyNA(w)
  if (has_na) {
    if (!na.rm) return(NULL)
    keep <- !is.na(x) & !is.na(y) & !is.na(w)
    x <- x[keep]
    y <- y[keep]
    w <- w[keep]
  }
  if (length(x) == 0L)
    rlang::abort("No non-NA observations after removing NAs.", call = call)
  list(x = x, y = y, w = w)
}

# ---- Exported functions ------------------------------------------------------

#' Weighted variance
#'
#' Computes the population (importance-weight) weighted variance
#' `sum(w_hat * (x - mu)^2)` where `w_hat = w / sum(w)` and
#' `mu = weighted.mean(x, w)`.
#'
#' @param x A numeric vector of observations.
#' @param w A numeric vector of non-negative weights (need not sum to 1).
#' @param na.rm Logical. If `TRUE`, paired `NA`s in `x` and `w` are removed
#'   before computation. If `FALSE` (default) and any `NA` is present,
#'   `NA` is returned.
#' @return A single numeric value.
#' @export
weighted_var <- function(x, w, na.rm = FALSE) {
  .check_wx(x, w)
  d <- .drop_na_wx(x, w, na.rm)
  if (is.null(d)) return(NA_real_)
  weighted_var_cpp(d$x, d$w)
}

#' Weighted standard deviation
#'
#' Square root of [weighted_var()].
#'
#' @inheritParams weighted_var
#' @return A single numeric value.
#' @export
weighted_sd <- function(x, w, na.rm = FALSE) {
  sqrt(weighted_var(x, w, na.rm = na.rm))
}

#' Weighted covariance
#'
#' Computes the population (importance-weight) weighted covariance
#' `sum(w_hat * (x - mu_x) * (y - mu_y))`.
#'
#' @param x A numeric vector of observations.
#' @param y A numeric vector of observations (same length as `x`).
#' @param w A numeric vector of non-negative weights (need not sum to 1).
#' @param na.rm Logical. If `TRUE`, observations with any `NA` in `x`, `y`,
#'   or `w` are removed before computation.
#' @return A single numeric value.
#' @export
weighted_cov <- function(x, y, w, na.rm = FALSE) {
  .check_wxy(x, y, w)
  d <- .drop_na_wxy(x, y, w, na.rm)
  if (is.null(d)) return(NA_real_)
  weighted_cov_cpp(d$x, d$y, d$w)
}

#' Weighted correlation
#'
#' Weighted correlation coefficient, computed as
#' `weighted_cov(x, y, w) / (weighted_sd(x, w) * weighted_sd(y, w))`.
#'
#' @inheritParams weighted_cov
#' @return A single numeric value in \[-1, 1\].
#' @export
weighted_cor <- function(x, y, w, na.rm = FALSE) {
  .check_wxy(x, y, w)
  d <- .drop_na_wxy(x, y, w, na.rm)
  if (is.null(d)) return(NA_real_)
  cov_xy <- weighted_cov_cpp(d$x, d$y, d$w)
  sd_x   <- sqrt(weighted_var_cpp(d$x, d$w))
  sd_y   <- sqrt(weighted_var_cpp(d$y, d$w))
  if (sd_x == 0 || sd_y == 0) return(NaN)
  cov_xy / (sd_x * sd_y)
}

#' Weighted quantiles
#'
#' Computes weighted quantiles using a type-7 analog: observations are placed
#' at the midpoints of their weight intervals in the cumulative weight
#' distribution, rescaled to \[0, 1\], and quantiles are obtained by linear
#' interpolation. For equal weights this matches `quantile(x, type = 7)`.
#'
#' @param x A numeric vector of observations.
#' @param w A numeric vector of non-negative weights (need not sum to 1).
#' @param probs A numeric vector of probabilities in \[0, 1\].
#' @param na.rm Logical. If `TRUE`, paired `NA`s in `x` and `w` are removed.
#' @return A numeric vector of the same length as `probs`.
#' @export
weighted_quantile <- function(x, w, probs = seq(0, 1, 0.25), na.rm = FALSE) {
  .check_wx(x, w)
  if (!is.numeric(probs) || any(probs < 0 | probs > 1, na.rm = TRUE))
    rlang::abort("`probs` must be a numeric vector with values in [0, 1].")
  d <- .drop_na_wx(x, w, na.rm)
  if (is.null(d)) return(rep(NA_real_, length(probs)))
  weighted_quantile_cpp(d$x, d$w, probs)
}

#' Weighted median
#'
#' Convenience wrapper: `weighted_quantile(x, w, 0.5, na.rm)`.
#'
#' @inheritParams weighted_var
#' @return A single numeric value.
#' @export
weighted_median <- function(x, w, na.rm = FALSE) {
  weighted_quantile(x, w, probs = 0.5, na.rm = na.rm)
}

#' Weighted median absolute deviation
#'
#' Computes the weighted MAD as the weighted median of `|x - median(x, w)|`,
#' scaled by `constant` (default 1.4826 for consistency with `stats::mad()`
#' under normality).
#'
#' @inheritParams weighted_var
#' @param constant Scale factor (default `1.4826`).
#' @return A single numeric value.
#' @export
weighted_mad <- function(x, w, na.rm = FALSE, constant = 1.4826) {
  .check_wx(x, w)
  d <- .drop_na_wx(x, w, na.rm)
  if (is.null(d)) return(NA_real_)
  med <- weighted_quantile_cpp(d$x, d$w, 0.5)
  constant * weighted_quantile_cpp(abs(d$x - med), d$w, 0.5)
}

#' Weighted skewness
#'
#' Computes the population weighted skewness (Fisher's g1):
#' `sum(w_hat * ((x - mu) / sigma)^3)`.
#'
#' @inheritParams weighted_var
#' @return A single numeric value.
#' @export
weighted_skewness <- function(x, w, na.rm = FALSE) {
  .check_wx(x, w)
  d <- .drop_na_wx(x, w, na.rm)
  if (is.null(d)) return(NA_real_)
  weighted_skewness_cpp(d$x, d$w)
}

#' Weighted excess kurtosis
#'
#' Computes the population weighted excess kurtosis:
#' `sum(w_hat * ((x - mu) / sigma)^4) - 3`.
#'
#' @inheritParams weighted_var
#' @return A single numeric value.
#' @export
weighted_kurtosis <- function(x, w, na.rm = FALSE) {
  .check_wx(x, w)
  d <- .drop_na_wx(x, w, na.rm)
  if (is.null(d)) return(NA_real_)
  weighted_kurtosis_cpp(d$x, d$w)
}
