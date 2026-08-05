#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

// Build the transition array for the augmented HSMM states.
// The R array has dimensions time x from x to.
// [[Rcpp::export]]
arma::cube compute_gamma(
    const arma::mat& omega,
    const arma::cube& omega_time,
    const arma::mat& p_array,
    const int K,
    const int M,
    const bool use_time_varying = false) {

  if (K < 1 || M < 1) Rcpp::stop("K and M must be positive.");

  const int Tm1 = p_array.n_rows;
  const int S = K * M;
  if (p_array.n_cols != static_cast<unsigned int>(S)) {
    Rcpp::stop("p_array must have K*M columns.");
  }
  if (!use_time_varying &&
      (omega.n_rows != static_cast<unsigned int>(K) ||
       omega.n_cols != static_cast<unsigned int>(K))) {
    Rcpp::stop("omega must be a K by K matrix.");
  }
  if (use_time_varying &&
      (omega_time.n_rows != static_cast<unsigned int>(Tm1) ||
       omega_time.n_cols != static_cast<unsigned int>(K) ||
       omega_time.n_slices != static_cast<unsigned int>(K))) {
    Rcpp::stop("omega_time must have dimensions time x K x K.");
  }

  arma::cube Gamma(Tm1, S, S, arma::fill::zeros);

  for (int t = 0; t < Tm1; ++t) {
    for (int k = 0; k < K; ++k) {
      for (int m = 0; m < M; ++m) {
        const int from = k * M + m;

        if (K == 1) {
          const int to = (m < M - 1) ? from + 1 : from;
          Gamma(t, from, to) = 1.0;
          continue;
        }

        double p = p_array(t, from);
        p = std::max(0.0, std::min(1.0, p));

        double weight_sum = 0.0;
        for (int i = 0; i < K; ++i) {
          if (i == k) continue;
          const double w = use_time_varying ?
            omega_time(t, k, i) : omega(k, i);
          if (std::isfinite(w) && w > 0.0) weight_sum += w;
        }

        const int stay_to = (m < M - 1) ? from + 1 : from;
        if (weight_sum <= 0.0) {
          Gamma(t, from, stay_to) = 1.0;
          continue;
        }

        Gamma(t, from, stay_to) = 1.0 - p;
        for (int i = 0; i < K; ++i) {
          if (i == k) continue;
          const double w = use_time_varying ?
            omega_time(t, k, i) : omega(k, i);
          if (std::isfinite(w) && w > 0.0) {
            Gamma(t, from, i * M) = p * w / weight_sum;
          }
        }
      }
    }
  }

  return Gamma;
}
