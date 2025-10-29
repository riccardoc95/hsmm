#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube compute_gamma(
    const arma::mat& omega,
    const arma::cube& omega_time,
    const arma::mat& p_array,
    const int K,
    const int M,
    const bool use_time_varying = false
) {
  const int n = p_array.n_rows;
  const int S = K * M;
  const int Tm1 = n - 1;

  arma::cube Gamma(n, S, S, arma::fill::zeros);
  arma::ivec dmap = arma::repelem(arma::regspace<arma::ivec>(0, K-1), M, 1);

  // parallelize on t (time)
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#pragma omp parallel for num_threads(nthreads)
#endif
  for (int t = 0; t < n; ++t) {
    for (int k = 0; k < K; ++k) {
      for (int m = 0; m < M-1; ++m) {
        int from_idx = k*M + m;
        int to_idx   = k*M + m + 1;
        Gamma(t, from_idx, to_idx) = 1.0 - p_array(t, from_idx);
      }

      int last_sub = k*M + (M - 1);
      for (int i = 0; i < K; ++i) {
        if (i == k) continue;
        int first_i = i*M;

        double omega_ki;
        if (use_time_varying) {
          int tt = std::min(t, Tm1 - 1);
          omega_ki = omega_time(tt, k, i);
        } else {
          omega_ki = omega(k, i);
        }

        Gamma(t, last_sub, first_i) = p_array(t, last_sub) * omega_ki;
      }
    }
  }

  return Gamma;
}
