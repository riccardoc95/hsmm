library(hsmm)

sim <- simulate_hsmm(
  n = 150,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 2,
  max_dwell = 10,
  seed = 123
)

fit <- fit_hsmm(
  sim$data,
  n_states = 3,
  model_type = "gaussian",
  max_dwell = 10,
  max_iter = 8,
  verbose = FALSE,
  seed = 123
)

stopifnot(
  is.finite(fit$final_loglik),
  nrow(fit$posterior) == nrow(sim$data),
  max(abs(rowSums(fit$posterior) - 1)) < 1e-7
)

print(fit)
