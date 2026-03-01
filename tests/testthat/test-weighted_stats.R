library(wstats)

# ---- helpers -----------------------------------------------------------------

n  <- 20L
set.seed(42)
x  <- rnorm(n)
y  <- rnorm(n)
w  <- rexp(n)          # unnormalised importance weights
w1 <- rep(1.0, n)      # equal weights

# Population variance / sd / cov from first principles
pop_var <- function(x, w = rep(1, length(x))) {
  wh <- w / sum(w)
  mu <- sum(wh * x)
  sum(wh * (x - mu)^2)
}
pop_cov <- function(x, y, w = rep(1, length(x))) {
  wh <- w / sum(w)
  mx <- sum(wh * x)
  my <- sum(wh * y)
  sum(wh * (x - mx) * (y - my))
}

# ---- weighted_var / weighted_sd ---------------------------------------------

test_that("weighted_var matches population variance formula", {
  expect_equal(weighted_var(x, w), pop_var(x, w))
})

test_that("weighted_var with equal weights equals population variance", {
  expect_equal(weighted_var(x, w1), pop_var(x))         # = var(x)*(n-1)/n
})

test_that("weighted_sd^2 == weighted_var", {
  expect_equal(weighted_sd(x, w)^2, weighted_var(x, w))
  expect_equal(weighted_sd(x, w1)^2, weighted_var(x, w1))
})

# ---- weighted_cov / weighted_cor --------------------------------------------

test_that("weighted_cov matches population covariance formula", {
  expect_equal(weighted_cov(x, y, w), pop_cov(x, y, w))
})

test_that("weighted_cov with equal weights equals population covariance", {
  expect_equal(weighted_cov(x, y, w1), pop_cov(x, y))   # = cov(x,y)*(n-1)/n
})

test_that("weighted_cor with equal weights equals cor(x, y)", {
  expect_equal(weighted_cor(x, y, w1), cor(x, y))
})

test_that("weighted_cor is in [-1, 1]", {
  r <- weighted_cor(x, y, w)
  expect_gte(r, -1)
  expect_lte(r,  1)
})

test_that("weighted_cor(x, x, w) == 1", {
  expect_equal(weighted_cor(x, x, w), 1)
})

# ---- weighted_quantile / weighted_median ------------------------------------

test_that("weighted_quantile with equal weights matches quantile(type=7)", {
  probs <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  expect_equal(
    weighted_quantile(x, w1, probs),
    unname(quantile(x, probs, type = 7))
  )
})

test_that("weighted_quantile(x, w, 0.5) == weighted_median(x, w)", {
  expect_equal(weighted_quantile(x, w, 0.5), weighted_median(x, w))
  expect_equal(weighted_quantile(x, w1, 0.5), weighted_median(x, w1))
})

test_that("weighted_quantile edge cases: p=0 gives min, p=1 gives max", {
  expect_equal(weighted_quantile(x, w, 0), min(x))
  expect_equal(weighted_quantile(x, w, 1), max(x))
  expect_equal(weighted_quantile(x, w1, 0), min(x))
  expect_equal(weighted_quantile(x, w1, 1), max(x))
})

test_that("weighted_quantile returns vector of correct length", {
  probs <- c(0.25, 0.5, 0.75)
  expect_length(weighted_quantile(x, w, probs), 3L)
})

# ---- weighted_mad -----------------------------------------------------------

test_that("weighted_mad with equal weights is close to mad(x)", {
  expect_equal(weighted_mad(x, w1), mad(x))
})

test_that("weighted_mad is non-negative", {
  expect_gte(weighted_mad(x, w), 0)
})

# ---- weighted_skewness / weighted_kurtosis ----------------------------------

test_that("weighted_skewness of symmetric data is near 0", {
  xs <- c(-3, -1, 0, 1, 3)
  ws <- rep(1, 5)
  expect_equal(weighted_skewness(xs, ws), 0, tolerance = 1e-10)
})

test_that("weighted_kurtosis of normal is near 0", {
  # For a large sample from N(0,1) with equal weights, excess kurtosis ~ 0
  set.seed(1)
  z <- rnorm(1e4)
  expect_equal(weighted_kurtosis(z, rep(1, 1e4)), 0, tolerance = 0.1)
})

test_that("weighted_kurtosis of uniform is near -1.2", {
  # Exact excess kurtosis of Uniform(-1,1) is -6/5 = -1.2
  set.seed(2)
  z <- runif(2e4, -1, 1)
  expect_equal(weighted_kurtosis(z, rep(1, 2e4)), -1.2, tolerance = 0.05)
})

# ---- NA handling ------------------------------------------------------------

test_that("functions return NA when NAs present and na.rm = FALSE", {
  xna <- c(x[1:5], NA, x[7:n])
  wna <- c(w[1:5], NA, w[7:n])

  expect_identical(weighted_var(xna, w),   NA_real_)
  expect_identical(weighted_var(x, wna),   NA_real_)
  expect_identical(weighted_sd(xna, w),    NA_real_)
  expect_identical(weighted_cov(xna, y, w), NA_real_)
  expect_identical(weighted_cov(x, y, wna), NA_real_)
  expect_identical(weighted_cor(xna, y, w), NA_real_)
  expect_equal(weighted_quantile(xna, w, 0.5), NA_real_)
  expect_identical(weighted_median(xna, w),  NA_real_)
  expect_identical(weighted_mad(xna, w),     NA_real_)
  expect_identical(weighted_skewness(xna, w), NA_real_)
  expect_identical(weighted_kurtosis(xna, w), NA_real_)
})

test_that("na.rm = TRUE gives same result as pre-filtering", {
  xna <- x; xna[3] <- NA
  wna <- w; wna[7] <- NA

  keep <- !is.na(xna) & !is.na(wna)
  xc   <- xna[keep]
  wc   <- wna[keep]

  expect_equal(weighted_var(xna, wna, na.rm = TRUE),  weighted_var(xc, wc))
  expect_equal(weighted_sd(xna, wna,  na.rm = TRUE),  weighted_sd(xc, wc))
  expect_equal(weighted_median(xna, wna, na.rm = TRUE), weighted_median(xc, wc))
  expect_equal(weighted_mad(xna, wna, na.rm = TRUE),  weighted_mad(xc, wc))
  expect_equal(weighted_skewness(xna, wna, na.rm = TRUE), weighted_skewness(xc, wc))
  expect_equal(weighted_kurtosis(xna, wna, na.rm = TRUE), weighted_kurtosis(xc, wc))
})

# ---- Input validation errors ------------------------------------------------

test_that("non-numeric x raises error", {
  expect_error(weighted_var("a", w),   "`x` must be a numeric vector.")
  expect_error(weighted_var(TRUE, w),  "`x` must be a numeric vector.")
})

test_that("non-numeric w raises error", {
  expect_error(weighted_var(x, "a"),   "`w` must be a numeric vector.")
})

test_that("length mismatch raises error with lengths in message", {
  expect_error(
    weighted_var(x, w[1:5]),
    regexp = "same length"
  )
})

test_that("negative weights raise error", {
  w_neg <- w; w_neg[1] <- -1
  expect_error(weighted_var(x, w_neg), "`w` must be non-negative.")
})

test_that("invalid probs raise error", {
  expect_error(weighted_quantile(x, w, probs = -0.1))
  expect_error(weighted_quantile(x, w, probs = 1.1))
  expect_error(weighted_quantile(x, w, probs = "a"))
})

test_that("weighted_cov/cor checks y length", {
  expect_error(weighted_cov(x, y[1:5], w), "same length")
  expect_error(weighted_cor(x, y[1:5], w), "same length")
})
