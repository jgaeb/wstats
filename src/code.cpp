#include <cpp11.hpp>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

using namespace cpp11;

// ---- internal helpers -------------------------------------------------------

static double wmean_impl(const std::vector<double>& x,
                         const std::vector<double>& w) {
  double sum_wx = 0.0, sum_w = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    sum_wx += w[i] * x[i];
    sum_w  += w[i];
  }
  return sum_wx / sum_w;
}

// Population (importance-weight) variance: sum(w_hat_i * (x_i - mu)^2)
// Input is assumed validated and NA-free by the R wrapper.
[[cpp11::register]]
double weighted_var_cpp(doubles x, doubles w) {
  int n = x.size();
  std::vector<double> xv(x.begin(), x.end());
  std::vector<double> wv(w.begin(), w.end());

  double V1 = 0.0;
  for (int i = 0; i < n; ++i) V1 += wv[i];

  double mu = wmean_impl(xv, wv);
  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    double d = xv[i] - mu;
    ss += wv[i] * d * d;
  }
  return ss / V1;
}

// Population (importance-weight) covariance: sum(w_hat_i * (x_i-mu_x) * (y_i-mu_y))
// Input is assumed validated and NA-free by the R wrapper.
[[cpp11::register]]
double weighted_cov_cpp(doubles x, doubles y, doubles w) {
  int n = x.size();
  std::vector<double> xv(x.begin(), x.end());
  std::vector<double> yv(y.begin(), y.end());
  std::vector<double> wv(w.begin(), w.end());

  double V1 = 0.0;
  for (int i = 0; i < n; ++i) V1 += wv[i];

  double mu_x = wmean_impl(xv, wv);
  double mu_y = wmean_impl(yv, wv);

  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    ss += wv[i] * (xv[i] - mu_x) * (yv[i] - mu_y);
  }
  return ss / V1;
}

// Weighted quantiles: type-7 analog.
//
// Sort (x, w) by x and normalise weights. Place each observation at the
// midpoint of its weight interval: M[i] = cumsum(w)[i] - w[i]/2. Rescale
// midpoints to [0, 1]: h[i] = (M[i] - M[0]) / (M[n-1] - M[0]). Linearly
// interpolate x at quantile p using h. For equal weights h[i] = i/(n-1),
// which matches R's quantile(type = 7). Edge cases: p <= h[0] -> x[0],
// p >= h[n-1] -> x[n-1].
//
// Input is assumed validated and NA-free by the R wrapper.
[[cpp11::register]]
doubles weighted_quantile_cpp(doubles x, doubles w, doubles probs) {
  int n = x.size();
  int np = probs.size();

  // Sort indices by x value
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&x](int a, int b) {
    return x[a] < x[b];
  });

  // Build sorted x and normalised weights
  std::vector<double> xs(n), ws(n);
  double sum_w = 0.0;
  for (int i = 0; i < n; ++i) {
    xs[i] = x[idx[i]];
    ws[i] = w[idx[i]];
    sum_w += ws[i];
  }
  for (int i = 0; i < n; ++i) ws[i] /= sum_w;

  // M[i] = midpoint of weight interval; then rescale to [0,1]
  std::vector<double> h(n);
  double cum = 0.0;
  for (int i = 0; i < n; ++i) {
    cum += ws[i];
    h[i] = cum - ws[i] / 2.0;
  }
  double h0 = h[0], hn = h[n - 1];
  double span = hn - h0;
  for (int i = 0; i < n; ++i) {
    h[i] = (span > 0.0) ? (h[i] - h0) / span : 0.0;
  }

  writable::doubles result(np);
  for (int j = 0; j < np; ++j) {
    double p = probs[j];
    if (n == 1 || p <= h[0]) {
      result[j] = xs[0];
    } else if (p >= h[n - 1]) {
      result[j] = xs[n - 1];
    } else {
      int lo = 0, hi = n - 1;
      while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (h[mid] <= p) lo = mid; else hi = mid;
      }
      double t = (p - h[lo]) / (h[hi] - h[lo]);
      result[j] = xs[lo] + t * (xs[hi] - xs[lo]);
    }
  }

  return result;
}

// Weighted skewness: sum(w_hat_i * ((x_i - mu) / sigma)^3)  [Fisher's g1]
// Input is assumed validated and NA-free by the R wrapper.
[[cpp11::register]]
double weighted_skewness_cpp(doubles x, doubles w) {
  int n = x.size();
  std::vector<double> xv(x.begin(), x.end());
  std::vector<double> wv(w.begin(), w.end());

  double V1 = 0.0;
  for (int i = 0; i < n; ++i) V1 += wv[i];

  double mu = wmean_impl(xv, wv);

  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    double d = xv[i] - mu;
    ss += wv[i] * d * d;
  }
  double sigma = std::sqrt(ss / V1);

  if (sigma == 0.0) return R_NaN;

  double s3 = 0.0;
  for (int i = 0; i < n; ++i) {
    double z = (xv[i] - mu) / sigma;
    s3 += wv[i] * z * z * z;
  }
  return s3 / V1;
}

// Weighted excess kurtosis: sum(w_hat_i * ((x_i - mu) / sigma)^4) - 3
// Input is assumed validated and NA-free by the R wrapper.
[[cpp11::register]]
double weighted_kurtosis_cpp(doubles x, doubles w) {
  int n = x.size();
  std::vector<double> xv(x.begin(), x.end());
  std::vector<double> wv(w.begin(), w.end());

  double V1 = 0.0;
  for (int i = 0; i < n; ++i) V1 += wv[i];

  double mu = wmean_impl(xv, wv);

  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    double d = xv[i] - mu;
    ss += wv[i] * d * d;
  }
  double sigma = std::sqrt(ss / V1);

  if (sigma == 0.0) return R_NaN;

  double s4 = 0.0;
  for (int i = 0; i < n; ++i) {
    double z = (xv[i] - mu) / sigma;
    s4 += wv[i] * z * z * z * z;
  }
  return s4 / V1 - 3.0;
}
