test_that("simulation and gaussian fit run end to end", {
  sim <- simulate_hsmm(
    n = 120, n_states = 3, n_dim = 2,
    model_type = "gaussian", max_dwell = 10, seed = 4
  )
  expect_s3_class(sim, "hsmm_simulation")
  expect_equal(dim(sim$data), c(120, 2))
  expect_true(all(sim$state %in% 1:3))

  fit <- fit_hsmm(
    sim$data, n_states = 3, model_type = "gaussian",
    max_dwell = 10, max_iter = 5, verbose = FALSE, seed = 4
  )
  expect_s3_class(fit, "hsmm_fit")
  expect_equal(dim(fit$posterior), c(120, 3))
  expect_equal(rowSums(fit$posterior), rep(1, 120), tolerance = 1e-8)
  expect_true(is.finite(fit$final_loglik))
  expect_equal(dim(fit$transition.model$gamma), c(119, 30, 30))

  new_sim <- simulate(fit, n = 30, seed = 5)
  expect_equal(nrow(new_sim$data), 30)
})

test_that("covariate models and HMM mode run", {
  set.seed(8)
  n <- 90
  x1 <- matrix(rnorm(n), n, 1)
  x2 <- cbind(sin(seq(0, 3 * pi, length.out = n)))
  sim <- simulate_hsmm(n, 3, covariates_q = x1,
                       covariates_omega = x2, seed = 8)

  fit <- fit_hsmm(
    sim$data, 3, covariates_q = x1, covariates_omega = x2,
    max_dwell = 7, max_iter = 3, verbose = FALSE, seed = 8
  )
  expect_equal(dim(fit$transition.model$omega_time), c(n - 1, 3, 3))
  expect_equal(ncol(fit$duration.model$beta.covariate), 1)

  hmm <- fit_hsmm(sim$data, 3, semi = FALSE, max_dwell = 20,
                  max_iter = 2, verbose = FALSE, seed = 8)
  expect_equal(hmm$input_params$max_dwell, 1)
})

test_that("all emission simulators return finite data", {
  families <- c("gaussian", "poisson", "exponential", "gamma",
                "beta", "student", "torus")#, "vargaussian")
  for (family in families) {
    sim <- simulate_hsmm(35, 2, model_type = family, n_dim = 2,
                         ar = 1, seed = 11)
    expect_true(all(is.finite(sim$data)), info = family)
  }
})

test_that("all emission families can be fitted for a few iterations", {
  families <- c("gaussian", "poisson", "exponential", "gamma",
                "beta", "student", "torus")#, "vargaussian")
  for (family in families) {
    sim <- simulate_hsmm(
      n = 55, n_states = 2, model_type = family,
      n_dim = 2, max_dwell = 5, ar = 1, seed = 20
    )
    fit <- fit_hsmm(
      sim$data, n_states = 2, model_type = family,
      max_dwell = 5, max_iter = 2, verbose = FALSE,
      ar = 1, seed = 20
    )
    expect_true(is.finite(fit$final_loglik), info = family)
    expect_equal(rowSums(fit$posterior), rep(1, 55),
                 tolerance = 1e-7, info = family)
  }
})

test_that("partial emission parameters use the remaining defaults", {
  custom_mu <- matrix(c(-6, 6), 2, 1)
  sim <- simulate_hsmm(
    n = 30, n_states = 2, model_type = "gaussian",
    emission = list(mu = custom_mu), seed = 9
  )
  expect_equal(sim$parameters$emission$mu, custom_mu)
  expect_length(sim$parameters$emission$sigma, 2)
})

test_that("simulation accepts transition-level covariates", {
  n <- 40
  x <- matrix(seq(-1, 1, length.out = n - 1), ncol = 1)
  sim <- simulate_hsmm(
    n, 3, covariates_q = x, covariates_omega = x, seed = 12
  )
  expect_equal(nrow(sim$covariates_q), n)
  expect_equal(nrow(sim$covariates_omega), n)
})

test_that("a supplied homogeneous omega is used directly", {
  omega <- matrix(c(
    0, 1, 0,
    1, 0, 0,
    0, 1, 0
  ), 3, 3, byrow = TRUE)
  sim <- simulate_hsmm(
    n = 50, n_states = 3, Pi = c(1, 0, 0), omega = omega,
    duration_intercept = rep(20, 3),
    duration_age_effect = rep(0, 3), seed = 2
  )
  expect_true(all(sim$state %in% c(1, 2)))
  expect_equal(sim$state[seq(1, 49, by = 2)], rep(1, 25))
  expect_equal(sim$state[seq(2, 50, by = 2)], rep(2, 25))
})


test_that("family is an alias in fitting and simulation", {
  sim <- simulate_hsmm(25, 2, family = "poisson", seed = 3)
  expect_true(all(sim$data >= 0))
  fit <- fit_hsmm(sim$data, 2, family = "poisson", max_dwell = 4,
                  max_iter = 1, verbose = FALSE, seed = 3)
  expect_equal(fit$input_params$model_type, "poisson")
})
