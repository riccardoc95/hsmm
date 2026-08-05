library(hsmm)
source("inst/tests/manual/test_helpers.R")


rates <- matrix(c(
  5.0, 5.0,
  1.0, 1.0,
  0.15, 0.15
), nrow = 3, byrow = TRUE)

sim <- simulate_hsmm(
  n = 900,
  n_states = 3,
  model_type = "exponential",
  n_dim = 2,
  max_dwell = 12,
  Pi = rep(1 / 3, 3),
  emission = list(rate = rates),
  seed = 106
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "exponential",
  max_dwell = 12,
  max_iter = 70,
  tol = 1e-4,
  verbose = FALSE,
  seed = 206
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Exponential, 3 states, 2 observed variables",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
