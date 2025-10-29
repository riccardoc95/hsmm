#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List backward_forward(const arma::mat& fit,
                                const arma::cube& Gamma,
                                const arma::vec& Pi_in,
                                const int K,
                                const int M) {


  const int n = fit.n_rows;
  const int S_expected = K * M;
  const int S_fit = fit.n_cols;

  if (n < 2) {
    stop("Need at least 2 observations (n >= 2).");
  }
  if (Gamma.n_rows != (n - 1)) {
    stop("Gamma must have n_obs-1 slices along first dim: dim(Gamma)[1] must be n-1.");
  }
  if (Gamma.n_cols != S_expected || Gamma.n_slices != S_expected) {
    stop("Gamma must have dimensions (n-1, K*M, K*M).");
  }
  if (Pi_in.n_elem != K) {
    stop("Pi must have length K.");
  }

  arma::vec Pi(S_expected, arma::fill::zeros);
  for (int k = 0; k < K; ++k) {
    Pi(k) = Pi_in(k); // first dwell slot gets Pi
  }

  arma::mat fit_expanded;
  if (S_fit == S_expected) {
    fit_expanded = fit;
  } else if (S_fit == K) {
    fit_expanded.set_size(n, S_expected);
    for (int k = 0; k < K; ++k) {
      for (int m = 0; m < M; ++m) {
        fit_expanded.col(k + m*K) = fit.col(k);
      }
    }
  } else {
    stop("fit must have either K or K*M columns.");
  }

  const int S = S_expected;
  arma::mat alpha(n, S, arma::fill::zeros);
  arma::mat Beta(n, S, arma::fill::zeros);
  arma::vec C(n, arma::fill::zeros);

  arma::rowvec psi(S, arma::fill::zeros);
  arma::rowvec psi_per_fit(S, arma::fill::zeros);


  psi = Pi.t();
  psi_per_fit = psi % fit_expanded.row(0);
  double c0 = arma::accu(psi_per_fit);
  if (c0 < DBL_EPSILON) c0 = DBL_EPSILON;
  C(0) = c0;
  alpha.row(0) = psi_per_fit / C(0);

  for (int i = 1; i < n; ++i) {

    arma::mat Gamma_mat(S, S, arma::fill::zeros);
    for (int h = 0; h < S; ++h) {
      for (int k2 = 0; k2 < S; ++k2) {
        Gamma_mat(h, k2) = Gamma(i-1, h, k2);
      }
    }

    psi = alpha.row(i-1) * Gamma_mat;
    psi_per_fit = psi % fit_expanded.row(i);

    double ci = arma::accu(psi_per_fit);
    if (ci < DBL_EPSILON) ci = DBL_EPSILON;
    C(i) = ci;

    alpha.row(i) = psi_per_fit / C(i);
  }

  double loglik = arma::sum(arma::log(C));

  Beta.row(n-1).ones();

  for (int i = n-2; i >= 0; --i) {

    arma::mat Gamma_mat(S, S, arma::fill::zeros);
    for (int h = 0; h < S; ++h) {
      for (int k2 = 0; k2 < S; ++k2) {
        Gamma_mat(h, k2) = Gamma(i, h, k2);
      }
    }

    arma::colvec tmp = (fit_expanded.row(i+1).t() % Beta.row(i+1).t());

    arma::colvec Bt = Gamma_mat * tmp;
    double denom = C(i+1);
    if (denom < DBL_EPSILON) denom = DBL_EPSILON;

    Beta.row(i) = (Bt/denom).t();
  }

  arma::mat post_pi_aug(n, S, arma::fill::zeros);
  arma::cube post_bi_pi_aug(n-1, S, S, arma::fill::zeros);

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#pragma omp parallel for num_threads(nthreads)
#endif
  for (int i = 0; i < n; ++i) {
    arma::rowvec num = alpha.row(i) % Beta.row(i);
    double z = arma::accu(num);
    if (z < DBL_EPSILON) z = DBL_EPSILON;
    post_pi_aug.row(i) = num / z;
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
  for (int i = 0; i < n-1; ++i) {

    double denom = C(i+1);
    if (denom < DBL_EPSILON) denom = DBL_EPSILON;

    for (int h = 0; h < S; ++h) {
      double ah = alpha(i, h);
      for (int k2 = 0; k2 < S; ++k2) {
        double trans_hk = Gamma(i, h, k2);
        double emit_k   = fit_expanded(i+1, k2);
        double beta_k   = Beta(i+1, k2);

        post_bi_pi_aug(i, h, k2) = (ah * trans_hk * emit_k * beta_k) / denom;
      }
    }

    double norm_i = 0.0;
    for (int h = 0; h < S; ++h) {
      for (int k2 = 0; k2 < S; ++k2) {
        norm_i += post_bi_pi_aug(i, h, k2);
      }
    }
    if (norm_i < DBL_EPSILON) norm_i = DBL_EPSILON;
    for (int h = 0; h < S; ++h) {
      for (int k2 = 0; k2 < S; ++k2) {
        post_bi_pi_aug(i, h, k2) /= norm_i;
      }
    }
  }

#ifdef _OPENMP
  int threads_used = nthreads;
#else
  int threads_used = 1;
#endif

  return Rcpp::List::create(
    _["l"]               = loglik,
    _["post.pi.aug"]     = post_pi_aug,
    _["post.bi.pi.aug"]  = post_bi_pi_aug,
    _["threads_used"]    = threads_used
  );
}
