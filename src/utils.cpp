#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::cube Gamma_f_cpp(const arma::mat& Omega, const arma::mat& p_array, int K, int M, arma::vec ld, arma::vec d) {
  
  int n = p_array.n_rows;
  arma::mat p = p_array;
  arma::cube Gamma(arma::sum(ld), arma::sum(ld), n, arma::fill::zeros);
  arma::uvec targetPos(K);
  targetPos(0) = 0;
  targetPos.subvec(1, K-1) = arma::find(arma::diff(d) == 1) + 1;
  arma::uvec pos; 
  int start, end;
  arma::mat app;
  
  for(int t = 0; t < n; t++) {
    for(int k = 0; k < K; k++) {
      pos = arma::find(d == k + 1);
      start = pos(0);
      end = pos(pos.n_elem-1);
      app = 1 - p.row(t).subvec(start, end);
      if(pos.n_elem >= 2) {
        if(pos.n_elem == 2) {
          Gamma.slice(t).submat(start, start + 1, end - 1, end) = app.row(0).subvec(0, app.n_elem-2);
        } else {
          Gamma.slice(t).submat(start, start + 1, end - 1, end) = arma::diagmat(app.row(0).subvec(0, app.n_elem-2));
        }
      }
      Gamma(end, end, t) = app(0, app.n_elem-1);
      
      for(int i = 0; i < K; i++) {
        if(i != k) {
          arma::mat app2 = p.row(t).subvec(start, end) * Omega(k, i);
          Gamma.slice(t).submat(start, targetPos(i), end, targetPos(i)) = app2.t();
        }
      }
    }
  }
  return Gamma;
}






// [[Rcpp::export]]
List lbackward_forward(arma::mat lfit, arma::cube Gamma, arma::vec Pi, arma::vec ld, int K, int M) {
  int n = lfit.n_rows;
  mat psi(n, sum(ld));
  mat Beta(n, sum(ld));
  mat alpha(n, sum(ld));
  vec C(n);
  mat post_pi_aug(n, sum(ld));
  double l;
  vec psi_per_fit(sum(ld));
  mat alphaBeta(n, sum(ld));
  cube post_bi_pi_aug(sum(ld), sum(ld), n - 1);
  
  // mat app  = repmat(Pi, M, 1);
  mat app  = Pi;
  psi.row(0) = app.t();
  psi_per_fit = vectorise(exp(log(psi.row(0)) + lfit.row(0)));
  C(0) = sum(psi_per_fit);
  alpha.row(0) = (exp(log(psi_per_fit) - log(C(0)))).t();
  
  for(int i = 1; i < n; i++) {
    psi.row(i) = alpha.row(i-1) * Gamma.slice(i - 1);
    psi_per_fit = vectorise(exp(log(psi.row(i)) + lfit.row(i)));
    C(i) = sum(psi_per_fit);
    alpha.row(i) = (exp(log(psi_per_fit) - log(C(i)))).t();
  }
  
  l = sum(log(C));
  
  for(int m = 0; m<sum(ld); m++){
    Beta(n-1, m) = 1;
  }
  for(int i = n - 2; i >= 0; i--) {
    Beta.row(i) = ((Gamma.slice(i) * (exp(lfit.row(i + 1) + log(Beta.row(i + 1)))).t()) / C(i+1)).t();
  }
  
  alphaBeta = exp(log(alpha) + log(Beta));
  for(int i = 0; i < n; i++){
    // post_pi_aug.row(i) = exp(log(alphaBeta.row(i)) - log(sum(alphaBeta.row(i))));
    post_pi_aug.row(i) = alphaBeta.row(i);
  }
  
  
  for(int i = 0; i < n - 1; i++) {
    for(int h = 0; h < sum(ld); h++) {
      for(int k = 0; k < sum(ld); k++) {
        post_bi_pi_aug(h, k, i) = exp(log(alpha(i,h)) + log(Gamma(h,k,i)) + lfit(i+1, k) + log(Beta(i+1, k)) - log(C(i+1)));
      }
    }
  }
  
  return List::create(Named("l") = l,
                      Named("post.pi.aug") = post_pi_aug,
                      Named("post.bi.pi.aug") = post_bi_pi_aug,
                      Named("alpha") = alpha,
                      Named("psi") = psi);
  
}
