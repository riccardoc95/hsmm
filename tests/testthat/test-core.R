test_that("augmented transition rows sum to one", {
  K <- 3
  M <- 4
  Tm1 <- 10
  omega <- matrix(1 / (K - 1), K, K)
  diag(omega) <- 0
  p <- matrix(runif(Tm1 * K * M, 0.05, 0.6), Tm1, K * M)

  G <- hsmm:::compute_gamma(
    omega, array(0, c(1, 1, 1)), p, K, M, FALSE
  )

  expect_equal(dim(G), c(Tm1, K * M, K * M))
  expect_equal(apply(G, c(1, 2), sum), matrix(1, Tm1, K * M),
               tolerance = 1e-12)
})

test_that("forward backward accepts state-level emissions", {
  K <- 2
  M <- 3
  n <- 12
  omega <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  p <- matrix(0.25, n - 1, K * M)
  G <- hsmm:::compute_gamma(
    omega, array(0, c(1, 1, 1)), p, K, M, FALSE
  )
  fit <- matrix(runif(n * K, 0.1, 1), n, K)
  ans <- hsmm:::backward.forward(fit, G, c(0.6, 0.4), K, M)

  expect_true(is.finite(ans$l))
  expect_equal(rowSums(ans$post.pi.aug), rep(1, n), tolerance = 1e-10)
  expect_equal(apply(ans$post.bi.pi.aug, 1, sum), rep(1, n - 1),
               tolerance = 1e-10)
})

test_that("conditional transition weights preserve the exit hazard", {
  K <- 3
  M <- 2
  p <- matrix(0.3, 1, K * M)
  omega <- matrix(c(
    4, 2, 6,
    3, 9, 1,
    2, 8, 5
  ), K, K, byrow = TRUE)

  G <- hsmm:::compute_gamma(
    omega, array(0, c(1, 1, 1)), p, K, M, FALSE
  )

  for (k in seq_len(K)) {
    from <- (k - 1) * M + 1
    same <- ((k - 1) * M + 1):(k * M)
    expect_equal(sum(G[1, from, -same]), 0.3, tolerance = 1e-12)
    expect_equal(sum(G[1, from, same]), 0.7, tolerance = 1e-12)
  }
})

test_that("initial probabilities start in the first dwell clone", {
  K <- 2
  M <- 3
  n <- 5
  omega <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  p <- matrix(0.2, n - 1, K * M)
  G <- hsmm:::compute_gamma(
    omega, array(0, c(1, 1, 1)), p, K, M, FALSE
  )
  fit <- matrix(1, n, K)
  ans <- hsmm:::backward.forward(fit, G, c(1, 0), K, M)

  expect_equal(ans$post.pi.aug[1, ], c(1, 0, 0, 0, 0, 0),
               tolerance = 1e-12)
})
