library(hsmm)
source("inst/tests/manual/test_helpers.R")


set.seed(111)
n <- 1000
x_q <- matrix(rnorm(n), ncol = 1)

sim <- simulate_hsmm(
  n = n,
  n_states = 2,
  model_type = "gaussian",
  n_dim = 1,
  max_dwell = 15,
  covariates_q = x_q,
  Pi = c(0.5, 0.5),
  duration_intercept = c(-2.7, -2.5),
  duration_age_effect = c(0.08, 0.08),
  duration_coef = matrix(c(1.2, -1.0), nrow = 2),
  emission = list(
    mu = matrix(c(-4, 4), nrow = 2),
    sigma = list(matrix(0.6, 1, 1), matrix(0.6, 1, 1))
  ),
  seed = 111
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 2,
  model_type = "gaussian",
  covariates_q = x_q,
  max_dwell = 15,
  max_iter = 80,
  tol = 1e-4,
  verbose = FALSE,
  seed = 211
)

test_result <- check_fit(
  sim,
  fit,
  minimum_accuracy = 0.95,
  test_name = "Gaussian with a duration covariate",
  make_plot = getOption("hsmm.make_plots", TRUE)
)

estimated_beta <- fit$duration.model$beta.covariate
hazard_variation <- max(apply(
  fit$duration.model$p.array,
  2,
  function(x) diff(range(x))
))

cat("Estimated duration coefficients:\n")
print(estimated_beta)
cat("Largest hazard variation over time:", round(hazard_variation, 3), "\n")

if (any(!is.finite(estimated_beta)) || hazard_variation < 0.01) {
  stop("The duration covariate does not appear to affect the fitted hazard")
}
