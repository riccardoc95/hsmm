library(hsmm)
source("inst/tests/manual/test_helpers.R")


sim <- simulate_hsmm(
  n = 900,
  n_states = 4,
  model_type = "gaussian",
  n_dim = 1,
  max_dwell = 12,
  Pi = rep(0.25, 4),
  duration_intercept = c(-2.6, -2.5, -2.4, -2.3),
  duration_age_effect = c(0.10, 0.10, 0.08, 0.08),
  emission = list(
    mu = matrix(c(-9, -3, 3, 9), nrow = 4),
    sigma = list(
      matrix(0.4, 1, 1),
      matrix(0.4, 1, 1),
      matrix(0.4, 1, 1),
      matrix(0.4, 1, 1)
    )
  ),
  seed = 103
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 4,
  model_type = "gaussian",
  max_dwell = 12,
  max_iter = 70,
  tol = 1e-4,
  verbose = FALSE,
  seed = 203
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gaussian, 4 states, 1 variable",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
