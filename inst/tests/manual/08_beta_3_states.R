library(hsmm)
source("inst/tests/manual/test_helpers.R")


means <- matrix(c(0.10, 0.50, 0.90), nrow = 3)
concentration <- 50

sim <- simulate_hsmm(
  n = 650,
  n_states = 3,
  model_type = "beta",
  n_dim = 1,
  max_dwell = 12,
  Pi = rep(1 / 3, 3),
  emission = list(
    shape1 = means * concentration,
    shape2 = (1 - means) * concentration
  ),
  seed = 108
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "beta",
  max_dwell = 12,
  max_iter = 70,
  tol = 1e-4,
  verbose = FALSE,
  seed = 208
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Beta, 3 states",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
