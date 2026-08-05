library(hsmm)
source("inst/tests/manual/test_helpers.R")


shape <- matrix(c(16, 16, 16), nrow = 3)
means <- matrix(c(0.8, 4, 12), nrow = 3)
rate <- shape / means

sim <- simulate_hsmm(
  n = 650,
  n_states = 3,
  model_type = "gamma",
  n_dim = 1,
  max_dwell = 12,
  Pi = rep(1 / 3, 3),
  emission = list(
    shape = shape,
    rate = rate
  ),
  seed = 107
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "gamma",
  max_dwell = 12,
  max_iter = 70,
  tol = 1e-4,
  verbose = FALSE,
  seed = 207
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gamma, 3 states",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
