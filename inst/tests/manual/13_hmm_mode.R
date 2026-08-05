library(hsmm)
source("inst/tests/manual/test_helpers.R")


sim <- simulate_hsmm(
  n = 600,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 1,
  semi = FALSE,
  max_dwell = 20,
  Pi = rep(1 / 3, 3),
  duration_intercept = c(-2.2, -2.0, -1.8),
  emission = list(
    mu = matrix(c(-5, 0, 5), nrow = 3),
    sigma = list(
      matrix(0.5, 1, 1),
      matrix(0.5, 1, 1),
      matrix(0.5, 1, 1)
    )
  ),
  seed = 113
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "gaussian",
  semi = FALSE,
  max_dwell = 20,
  max_iter = 60,
  tol = 1e-4,
  verbose = FALSE,
  seed = 213
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gaussian HMM mode, 3 states",
  make_plot = getOption("hsmm.make_plots", TRUE)
)

cat("max_dwell used by the fit:", fit$input_params$max_dwell, "\n")

if (fit$input_params$max_dwell != 1) {
  stop("semi = FALSE did not force max_dwell to 1")
}
