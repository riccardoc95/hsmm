library(hsmm)
source("inst/tests/manual/test_helpers.R")


sim <- simulate_hsmm(
  n = 600,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 1,
  max_dwell = 12,
  Pi = c(1 / 3, 1 / 3, 1 / 3),
  duration_intercept = c(-2.6, -2.4, -2.2),
  duration_age_effect = c(0.12, 0.10, 0.08),
  emission = list(
    mu = matrix(c(-6, 0, 6), nrow = 3),
    sigma = list(
      matrix(0.45, 1, 1),
      matrix(0.45, 1, 1),
      matrix(0.45, 1, 1)
    )
  ),
  seed = 102
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "gaussian",
  max_dwell = 12,
  max_iter = 60,
  tol = 1e-4,
  verbose = FALSE,
  seed = 202
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gaussian, 3 states, 1 variable",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
