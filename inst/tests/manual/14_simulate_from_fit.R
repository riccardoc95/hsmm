library(hsmm)
source("inst/tests/manual/test_helpers.R")


first_sim <- simulate_hsmm(
  n = 650,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 2,
  max_dwell = 12,
  emission = list(
    mu = matrix(c(
      -5, -5,
       0,  5,
       5, -5
    ), nrow = 3, byrow = TRUE),
    sigma = list(diag(2) * 0.5, diag(2) * 0.5, diag(2) * 0.5)
  ),
  seed = 114
)

first_fit <- fit_hsmm(
  data = first_sim$data,
  n_states = 3,
  model_type = "gaussian",
  max_dwell = 12,
  max_iter = 60,
  tol = 1e-4,
  verbose = FALSE,
  seed = 214
)

first_result <- check_fit(
  first_sim,
  first_fit,
  minimum_accuracy = 0.88,
  test_name = "Original fit before simulate(fit)",
  make_plot = FALSE
)

second_sim <- simulate(first_fit, n = 450, seed = 314)

second_fit <- fit_hsmm(
  data = second_sim$data,
  n_states = 3,
  model_type = "gaussian",
  max_dwell = 12,
  max_iter = 60,
  tol = 1e-4,
  verbose = FALSE,
  seed = 414
)

test_result <- check_fit(
  second_sim,
  second_fit,
  minimum_accuracy = 0.95,
  test_name = "Simulation from a fitted model",
  make_plot = getOption("hsmm.make_plots", TRUE)
)
