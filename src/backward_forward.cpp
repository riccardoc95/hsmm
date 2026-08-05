#include <RcppArmadillo.h>
#include <cfloat>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

static arma::mat gamma_at(const arma::cube& Gamma, int t, int S) {
  arma::mat out(S, S, arma::fill::zeros);
  for (int from = 0; from < S; ++from) {
    for (int to = 0; to < S; ++to) {
      out(from, to) = Gamma(t, from, to);
    }
  }
  return out;
}

// Scaled forward-backward algorithm for the augmented state space.
// [[Rcpp::export]]
Rcpp::List backward_forward(
    const arma::mat& fit,
    const arma::cube& Gamma,
    const arma::vec& Pi_in,
    const int K,
    const int M) {

  const int n = fit.n_rows;
  const int S = K * M;

  if (n < 2) Rcpp::stop("At least two observations are needed.");
  if (K < 1 || M < 1) Rcpp::stop("K and M must be positive.");
  if (Gamma.n_rows != static_cast<unsigned int>(n - 1) ||
      Gamma.n_cols != static_cast<unsigned int>(S) ||
      Gamma.n_slices != static_cast<unsigned int>(S)) {
    Rcpp::stop("Gamma must have dimensions (n-1) x (K*M) x (K*M).");
  }
  if (Pi_in.n_elem != static_cast<unsigned int>(K)) {
    Rcpp::stop("Pi must have length K.");
  }

  arma::mat fit_aug;
  if (fit.n_cols == static_cast<unsigned int>(S)) {
    fit_aug = fit;
  } else if (fit.n_cols == static_cast<unsigned int>(K)) {
    fit_aug.set_size(n, S);
    for (int k = 0; k < K; ++k) {
      for (int m = 0; m < M; ++m) {
        fit_aug.col(k * M + m) = fit.col(k);
      }
    }
  } else {
    Rcpp::stop("fit must have K or K*M columns.");
  }

  arma::rowvec Pi(S, arma::fill::zeros);
  for (int k = 0; k < K; ++k) {
    Pi(k * M) = Pi_in(k);
  }

  arma::mat alpha(n, S, arma::fill::zeros);
  arma::mat beta(n, S, arma::fill::zeros);
  arma::vec scale(n, arma::fill::zeros);

  arma::rowvec value = Pi % fit_aug.row(0);
  scale(0) = arma::accu(value);
  if (!std::isfinite(scale(0)) || scale(0) < DBL_MIN) scale(0) = DBL_MIN;
  alpha.row(0) = value / scale(0);

  for (int t = 1; t < n; ++t) {
    arma::mat G = gamma_at(Gamma, t - 1, S);
    value = (alpha.row(t - 1) * G) % fit_aug.row(t);
    scale(t) = arma::accu(value);
    if (!std::isfinite(scale(t)) || scale(t) < DBL_MIN) scale(t) = DBL_MIN;
    alpha.row(t) = value / scale(t);
  }

  const double loglik = arma::sum(arma::log(scale));
  beta.row(n - 1).ones();

  for (int t = n - 2; t >= 0; --t) {
    arma::mat G = gamma_at(Gamma, t, S);
    arma::colvec next = fit_aug.row(t + 1).t() % beta.row(t + 1).t();
    beta.row(t) = (G * next / scale(t + 1)).t();
  }

  arma::mat post_pi_aug(n, S, arma::fill::zeros);
  for (int t = 0; t < n; ++t) {
    arma::rowvec num = alpha.row(t) % beta.row(t);
    double total = arma::accu(num);
    if (!std::isfinite(total) || total < DBL_MIN) total = DBL_MIN;
    post_pi_aug.row(t) = num / total;
  }

  arma::cube post_bi_pi_aug(n - 1, S, S, arma::fill::zeros);
  for (int t = 0; t < n - 1; ++t) {
    double total = 0.0;
    for (int from = 0; from < S; ++from) {
      for (int to = 0; to < S; ++to) {
        const double x = alpha(t, from) * Gamma(t, from, to) *
          fit_aug(t + 1, to) * beta(t + 1, to) / scale(t + 1);
        post_bi_pi_aug(t, from, to) = x;
        total += x;
      }
    }
    if (!std::isfinite(total) || total < DBL_MIN) total = DBL_MIN;
    for (int from = 0; from < S; ++from) {
      for (int to = 0; to < S; ++to) {
        post_bi_pi_aug(t, from, to) /= total;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_["l"] = loglik,
    Rcpp::_["post.pi.aug"] = post_pi_aug,
    Rcpp::_["post.bi.pi.aug"] = post_bi_pi_aug
  );
}
