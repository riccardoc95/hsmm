library(hsmm)
source("inst/tests/manual/test_helpers.R")


sim <- simulate_hsmm(
  n = 700,
  n_states = 3,
  model_type = "torus",
  n_dim = 1,
  max_dwell = 12,
  Pi = rep(1 / 3, 3),
  emission = list(
    mu = matrix(c(0, 2 * pi / 3, 4 * pi / 3), nrow = 3),
    kappa = matrix(c(25, 25, 25), nrow = 3)
  ),
  seed = 110
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "torus",
  max_dwell = 12,
  max_iter = 70,
  tol = 1e-4,
  verbose = FALSE,
  seed = 210
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Torus, 3 states",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
