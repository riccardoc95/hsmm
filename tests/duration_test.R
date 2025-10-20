test_duration_models <- function() {
  set.seed(123)
  n_obs <- 100
  n_states <- 3
  max_dwell <- 5
  data <- matrix(rnorm(n_obs * 3), ncol = 3)
  
  simulate_post_bi_pi_aug <- function(n_obs, n_states, dwell_lengths) {
    total_dwell <- sum(dwell_lengths)
    array(runif((n_obs - 1) * total_dwell * total_dwell, 0, 1), 
          dim = c(n_obs - 1, total_dwell, total_dwell))
  }
  
  check_test <- function(model, params, post) {
    model$compute_p.array(params)
    model$compute_beta(params, post)
    
    p_array_ok <- all(dim(model$p.array) == c(params$n_obs - 1, sum(params$dwell_lengths)))
    beta_intercepts_ok <- length(model$beta.intercepts) == params$n_states
    beta_time_effects_ok <- length(model$beta.time_effects) == params$n_states
    beta_covar_ok <- is.null(model$beta.covariate) || ncol(model$beta.covariate) == ifelse(is.null(params$covariates_q), 0, ncol(params$covariates_q))
    
    return(p_array_ok && beta_intercepts_ok && beta_time_effects_ok && beta_covar_ok)
  }
  
  # Caso 1: Geometric no covariates, semi = FALSE
  params1 <- list(
    n_states = n_states,
    n_obs = n_obs,
    max_dwell = max_dwell,
    dwell_lengths = rep(max_dwell, n_states),
    covariates_q = NULL,
    semi = FALSE,
    state_indices = rep(1:n_states, each = max_dwell)
  )
  geom1 <- DurationModel.Geometric$new(params1)
  post1 <- simulate_post_bi_pi_aug(n_obs, n_states, params1$dwell_lengths)
  if (check_test(geom1, params1, post1)) {
    cat("Test 1 compute_p.array e compute_beta superati\n")
  } else {
    cat("Errore Test 1\n")
  }
  
  # Caso 2: Geometric no covariates, semi = TRUE
  params2 <- params1
  params2$semi <- TRUE
  geom2 <- DurationModel.Geometric$new(params2)
  post2 <- simulate_post_bi_pi_aug(n_obs, n_states, params2$dwell_lengths)
  if (check_test(geom2, params2, post2)) {
    cat("Test 2 compute_p.array e compute_beta superati\n")
  } else {
    cat("Errore Test 2\n")
  }
  
  # Caso 3: Gompertz covariates, semi = TRUE
  cov_q <- matrix(rnorm((n_obs - 1) * 2), ncol = 2)
  params3 <- list(
    n_states = n_states,
    n_obs = n_obs,
    max_dwell = max_dwell,
    dwell_lengths = rep(max_dwell, n_states),
    covariates_q = cov_q,
    semi = TRUE,
    state_indices = rep(1:n_states, each = max_dwell)
  )
  gomp <- DurationModel.Gompertz$new(params3)
  post3 <- simulate_post_bi_pi_aug(n_obs, n_states, params3$dwell_lengths)
  if (check_test(gomp, params3, post3)) {
    cat("Test 3 compute_p.array e compute_beta superati\n")
  } else {
    cat("Errore Test 3\n")
  }
  
  # Caso 4: Gompertz covariates, semi = FALSE
  params4 <- params3
  params4$semi <- FALSE
  gomp2 <- DurationModel.Gompertz$new(params4)
  post4 <- simulate_post_bi_pi_aug(n_obs, n_states, params4$dwell_lengths)
  if (check_test(gomp2, params4, post4)) {
    cat("Test 4 compute_p.array e compute_beta superati\n")
  } else {
    cat("Errore Test 4\n")
  }
}

test_duration_models()
