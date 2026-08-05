library(hsmm)
source("inst/tests/manual/test_helpers.R")


sim <- simulate_hsmm(
  n = 650,
  n_states = 3,
  model_type = "poisson",
  n_dim = 1,
  max_dwell = 12,
  Pi = rep(1 / 3, 3),
  emission = list(
    lambda = matrix(c(2, 10, 25), nrow = 3)
  ),
  seed = 105
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "poisson",
  max_dwell = 12,
  max_iter = 60,
  tol = 1e-4,
  verbose = FALSE,
  seed = 205
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Poisson, 3 states",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
