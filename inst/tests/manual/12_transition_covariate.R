library(hsmm)
source("inst/tests/manual/test_helpers.R")


n <- 1400
x_omega <- matrix(sin(seq(0, 8 * pi, length.out = n)), ncol = 1)

transition_coef <- list(
  matrix(c(0, 2.0), nrow = 2),
  matrix(c(0, -2.0), nrow = 2),
  matrix(c(0, 1.5), nrow = 2)
)

sim <- simulate_hsmm(
  n = n,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 1,
  max_dwell = 12,
  covariates_omega = x_omega,
  transition_coef = transition_coef,
  Pi = rep(1 / 3, 3),
  duration_intercept = c(-2.2, -2.2, -2.2),
  duration_age_effect = c(0.08, 0.08, 0.08),
  emission = list(
    mu = matrix(c(-6, 0, 6), nrow = 3),
    sigma = list(
      matrix(0.5, 1, 1),
      matrix(0.5, 1, 1),
      matrix(0.5, 1, 1)
    )
  ),
  seed = 112
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "gaussian",
  covariates_omega = x_omega,
  max_dwell = 12,
  max_iter = 80,
  tol = 1e-4,
  verbose = FALSE,
  seed = 212
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gaussian with a transition covariate",
  make_plot = getOption("hsmm.make_plots", TRUE)
)

all_coefficients <- unlist(fit$transition.model$coefficients)
omega_variation <- max(apply(
  fit$transition.model$omega_time,
  c(2, 3),
  function(x) diff(range(x))
))

cat("Estimated transition coefficients:\n")
print(fit$transition.model$coefficients)
cat("Largest transition probability variation:", round(omega_variation, 3), "\n")

if (any(!is.finite(all_coefficients)) || omega_variation < 0.01) {
  stop("The transition covariate does not appear to affect the fitted transitions")
}
