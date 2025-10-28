#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List backward_forward(const arma::mat& fit,      // n x S? (or n x K, but we expand)
                                const arma::cube& Gamma,   // (n-1) x S x S  <-- NOTE: R stores [t, h, k]
                                const arma::vec& Pi_in,    // length K
                                const int K,
                                const int M) {

  // -------------------------------------------------
  // Dimensions and sanity checks
  // -------------------------------------------------
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
    // Careful: user said Gamma is (n_obs-1, S, S),
    // but Armadillo cube is (n_rows, n_cols, n_slices),
    // so we interpret Gamma as (n-1, S, S) = (time, from, to).
    // That means n_rows = n-1, n_cols = S, n_slices = S.
    // We'll handle transpose access below.
    stop("Gamma must have dimensions (n-1, K*M, K*M).");
  }
  if (Pi_in.n_elem != K) {
    stop("Pi must have length K.");
  }

  // -------------------------------------------------
  // Expand Pi (length K) to Pi_expanded (length S=K*M)
  // Each macro-state k gets replicated M times
  // exactly like in R:
  // post.pi.init <- c(matrix(c(Pi, rep(0, (M-1)*K)), ncol=K, nrow=M, byrow=TRUE))
  // -------------------------------------------------
  arma::vec Pi(S_expected, arma::fill::zeros);
  for (int k = 0; k < K; ++k) {
    Pi(k) = Pi_in(k); // first dwell slot gets Pi
  }
  // NOTE: the rest (k + m*K for m>0) remain 0, consistent
  // with the HSMM augmentation

  // -------------------------------------------------
  // If fit has only K columns, expand it to S=K*M by
  // replicating each state's emission prob M times.
  // User said fit is already n x (K*M), but let's be safe.
  // -------------------------------------------------
  arma::mat fit_expanded;
  if (S_fit == S_expected) {
    fit_expanded = fit;                // already S states
  } else if (S_fit == K) {
    // replicate each column M times
    fit_expanded.set_size(n, S_expected);
    for (int k = 0; k < K; ++k) {
      for (int m = 0; m < M; ++m) {
        fit_expanded.col(k + m*K) = fit.col(k);
      }
    }
  } else {
    stop("fit must have either K or K*M columns.");
  }

  const int S = S_expected; // shorthand

  // -------------------------------------------------
  // We'll treat Gamma as time-indexed:
  // Gamma[t, h, k] in R  -> probability h->k from t to t+1
  //
  // In Armadillo, our cube is (n-1, S, S)
  //   => row = t, col = h, slice = k
  //
  // But Armadillo expects (n_rows, n_cols, n_slices).
  // We got (n-1, S, S).
  // We'll access Gamma_t(h,k) by Gamma(t, h, k).
  //
  // We'll build on-the-fly a arma::mat Gamma_t of shape SxS for speed.
  // -------------------------------------------------

  // -------------------------------------------------
  // Allocate forward/backward buffers
  // -------------------------------------------------
  arma::mat alpha(n, S, arma::fill::zeros);
  arma::mat Beta(n, S, arma::fill::zeros);
  arma::vec C(n, arma::fill::zeros);

  // temporary buffers
  arma::rowvec psi(S, arma::fill::zeros);
  arma::rowvec psi_per_fit(S, arma::fill::zeros);

  // -------------------------------------------------
  // Forward pass
  // -------------------------------------------------

  // t = 0 init
  psi = Pi.t(); // Pi is length S
  psi_per_fit = psi % fit_expanded.row(0);
  double c0 = arma::accu(psi_per_fit);
  if (c0 < DBL_EPSILON) c0 = DBL_EPSILON;
  C(0) = c0;
  alpha.row(0) = psi_per_fit / C(0);

  // t = 1..n-1
  for (int i = 1; i < n; ++i) {

    // Build Gamma_(i-1) matrix SxS from the cube
    arma::mat Gamma_mat(S, S, arma::fill::zeros);
    // Gamma_mat(h,k) = Gamma(i-1, h, k)
    // parallel fill
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

  // log-likelihood = sum log C
  double loglik = arma::sum(arma::log(C));

  // -------------------------------------------------
  // Backward pass
  // -------------------------------------------------

  // last row:
  // in R: Beta[n, ] <- 1
  // when using scaling C, strictly it's Beta[n,] <- 1
  // and at step i we divide by C[i+1]
  Beta.row(n-1).ones();

  // i from n-2 down to 0
  for (int i = n-2; i >= 0; --i) {

    arma::mat Gamma_mat(S, S, arma::fill::zeros);
    for (int h = 0; h < S; ++h) {
      for (int k2 = 0; k2 < S; ++k2) {
        Gamma_mat(h, k2) = Gamma(i, h, k2);
      }
    }

    // compute Beta[i,] = Gamma[i,,] %*% (fit[i+1,]*Beta[i+1,]) / C[i+1]
    arma::colvec tmp = (fit_expanded.row(i+1).t() % Beta.row(i+1).t());

    arma::colvec Bt = Gamma_mat * tmp;
    double denom = C(i+1);
    if (denom < DBL_EPSILON) denom = DBL_EPSILON;

    Beta.row(i) = (Bt/denom).t();
  }

  // -------------------------------------------------
  // post.pi.aug (posterior over states at each time)
  // post.bi.pi.aug (posterior over transitions)
  // -------------------------------------------------

  arma::mat post_pi_aug(n, S, arma::fill::zeros);
  arma::cube post_bi_pi_aug(n-1, S, S, arma::fill::zeros);
  // NOTE: shape matches R output [n-1, S, S]

  // univariate posterior
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

  // bivariate posterior
  // for each t, for each (h,k):
  // alpha[t,h] * Gamma[t,h,k] * fit[t+1,k] * Beta[t+1,k] / C[t+1]
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
  for (int i = 0; i < n-1; ++i) {

    double denom = C(i+1);
    if (denom < DBL_EPSILON) denom = DBL_EPSILON;

    for (int h = 0; h < S; ++h) {
      double ah = alpha(i, h);
      for (int k2 = 0; k2 < S; ++k2) {
        // Gamma(i, h, k2) from cube
        double trans_hk = Gamma(i, h, k2);
        double emit_k   = fit_expanded(i+1, k2);
        double beta_k   = Beta(i+1, k2);

        post_bi_pi_aug(i, h, k2) = (ah * trans_hk * emit_k * beta_k) / denom;
      }
    }

    // normalize slice i so rowsum over (h,k) = 1 (like R code effectively does with / rowSums(...) per-time)
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

  // -------------------------------------------------
  // threads_used
  // -------------------------------------------------
#ifdef _OPENMP
  int threads_used = nthreads;
#else
  int threads_used = 1;
#endif

  // -------------------------------------------------
  // Return
  // -------------------------------------------------
  return Rcpp::List::create(
    _["l"]               = loglik,
    _["post.pi.aug"]     = post_pi_aug,
    _["post.bi.pi.aug"]  = post_bi_pi_aug,
    _["threads_used"]    = threads_used
  );
}
