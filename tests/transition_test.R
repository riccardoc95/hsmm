library(testthat)
library(magrittr)
library(LaplacesDemon)
library(dplyr)
library(tidyr)
library(nnet)

# Funzione per generare dati fittizi per test parametri
generate_params <- function(n_states, n_obs, covariate_flag = FALSE, max_dwell = 5) {
  # Durate massime per ogni stato
  dwell_lengths <- rep(max_dwell, n_states)
  
  # Vettore d: indicizza le colonne di p.array secondo stato
  d <- rep(1:n_states, dwell_lengths)
  
  # Numero di colonne p.array = somma delle durate per ogni stato
  n_cols <- sum(dwell_lengths)
  
  # p.array: inizializzazione con zeri o valori casuali
  p.array <- array(0, dim = c(n_obs - 1, n_cols))
  
  # Per esempio, assegna valori random per test
  for (t in 1:(n_obs - 1)) {
    for (k in 1:n_states) {
      pos <- which(d == k)  # colonne corrispondenti allo stato k
      # assegna un vettore casuale come placeholder 
      # (qui puoi sostituire con valore reale di hazard/p_k)
      p.array[t, pos] <- runif(length(pos), min = 0, max = 1)
    }
  }
  
  params <- list(
    n_states = n_states,
    n_obs = n_obs,
    max_dwell = max_dwell,
    dwell_lengths = dwell_lengths,
    state_indices = rep(1:n_states, each = max_dwell),  # o altro a seconda del caso
    p.array = p.array,
    covariates = if (covariate_flag) matrix(rnorm((n_obs - 1) * 2), ncol = 2) else NULL
  )
  
  return(params)
}


test_transition_models_detailed <- function() {
  configs <- expand.grid(
    n_states = c(2, 3, 5),
    covariates = c(TRUE, FALSE)
  )
  
  for(row in seq_len(nrow(configs))) {
    n_states <- configs$n_states[row]
    covariate_flag <- configs$covariates[row]
    
    cat(sprintf("\nTesting n_states = %d, covariates = %s\n", n_states, covariate_flag))
    
    params <- generate_params(n_states, n_obs = 50, covariate_flag = covariate_flag)
    
    # TransitionModel.Homogeneous
    if (!covariates){
      cat("Test TransitionModel.Homogeneous initialize...\n")
      hom_model <- TransitionModel.Homogeneous$new(params, params$p.array)
      expect_equal(dim(hom_model$omega), c(n_states, n_states))
      expect_true(all(rowSums(hom_model$omega) == 1))
      
      cat("Test TransitionModel.Homogeneous updategamma...\n")
      hom_model$compute_gamma(params, params$p.array)
      expect_true(!is.null(hom_model$gamma))
      
      # Create dummy post.bi.pi.aug array with matching dimensions
      post_bi_pi_aug <- array(dim = c(params$n_obs-1, sum(params$dwell_lengths), sum(params$dwell_lengths)))
      #post_bi_pi_aug <- array(runif((params$n_obs-1)*n_states*n_states), dim = c(params$n_obs-1, n_states, n_states))
      post_bi_pi_aug[is.na(post_bi_pi_aug)] <- 1
  
      cat("Test TransitionModel.Homogeneous updateomega...\n")
      hom_model$compute_omega(params, post.bi.pi.aug = post_bi_pi_aug)
      expect_equal(dim(hom_model$omega), c(n_states, n_states))
      expect_true(all(rowSums(hom_model$omega) > 0)) # normalized rows
    } else{
      # TransitionModel.TimeVarying
      cat("Test TransitionModel.TimeVarying initialize...\n")
      time_model <- TransitionModel.TimeVarying$new(params, params$p.array)
      expect_equal(length(dim(time_model$omega)), 3)
      expect_equal(dim(time_model$omega)[2], n_states)
      expect_equal(dim(time_model$omega)[3], n_states)
      
      cat("Test TransitionModel.TimeVarying updategamma...\n")
      time_model$compute_gamma(params, params$p.array)
      expect_true(!is.null(time_model$gamma))
      
      cat("Test TransitionModel.TimeVarying updateomega...\n")
      time_model$compute_omega(params, covariates = params$covariates, post.bi.pi.aug = post_bi_pi_aug)
      expect_equal(length(dim(time_model$omega)), 3)
      expect_equal(dim(time_model$omega)[2], n_states)
      expect_equal(dim(time_model$omega)[3], n_states)
      for(t in 1:(length(dim(time_model$omega)) > 0 && dim(time_model$omega)[1])) {
        expect_true(abs(sum(time_model$omega[t, 1, ]) - 1) < 1e-6 ||
                      abs(sum(time_model$omega[t, 2, ]) - 1) < 1e-6)
      }
    }
    
    cat("All tests passed for configuration n_states =", n_states,
        "covariates =", covariate_flag, "\n")
  }
  cat("\nAll transition model tests completed successfully.\n")
}

# Per eseguire i test
test_transition_models_detailed()
