library(hsmm)
source("inst/tests/manual/test_helpers.R")

sim <- simulate_hsmm(
  n = 400,
  n_states = 2,
  model_type = "gaussian",
  n_dim = 1,
  max_dwell = 10,
  Pi = c(0.5, 0.5),
  duration_intercept = c(-2.4, -2.2),
  duration_age_effect = c(0.18, 0.14),
  emission = list(
    mu = matrix(c(-4, 4), nrow = 2),
    sigma = list(matrix(0.5, 1, 1), matrix(0.5, 1, 1))
  ),
  seed = 101
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 2,
  model_type = "gaussian",
  max_dwell = 10,
  max_iter = 50,
  tol = 1e-4,
  verbose = FALSE,
  seed = 201
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gaussian, 2 states, 1 variable",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
