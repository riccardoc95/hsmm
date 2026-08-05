library(hsmm)
source("inst/tests/manual/test_helpers.R")


mu <- matrix(c(
  -5, -5,
   0,  5,
   5, -5
), nrow = 3, byrow = TRUE)

sim <- simulate_hsmm(
  n = 650,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 2,
  max_dwell = 12,
  Pi = rep(1 / 3, 3),
  emission = list(
    mu = mu,
    sigma = list(diag(2) * 0.5, diag(2) * 0.5, diag(2) * 0.5)
  ),
  seed = 104
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "gaussian",
  max_dwell = 12,
  max_iter = 60,
  tol = 1e-4,
  verbose = FALSE,
  seed = 204
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gaussian, 3 states, 2 observed variables",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
