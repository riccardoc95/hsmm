library(hsmm)
source("inst/tests/manual/test_helpers.R")


sim <- simulate_hsmm(
  n = 700,
  n_states = 3,
  model_type = "student",
  n_dim = 1,
  max_dwell = 12,
  Pi = rep(1 / 3, 3),
  emission = list(
    mu = matrix(c(-6, 0, 6), nrow = 3),
    scale = matrix(c(0.6, 0.6, 0.6), nrow = 3),
    df = c(7, 7, 7)
  ),
  seed = 109
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "student",
  max_dwell = 12,
  max_iter = 70,
  tol = 1e-4,
  verbose = FALSE,
  seed = 209
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Student t, 3 states",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
