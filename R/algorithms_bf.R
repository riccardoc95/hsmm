# R wrapper around the C++ scaled forward-backward algorithm.
backward.forward <- function(fit, Gamma, Pi, K, M) {
  out <- backward_forward(fit, Gamma, Pi, K, M)
  list(
    l = out$l,
    post.pi.aug = out$post.pi.aug,
    post.bi.pi.aug = out$post.bi.pi.aug
  )
}
