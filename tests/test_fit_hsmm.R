# ============================================================================
# Script di Test per la funzione fit_hsmm
# ============================================================================
# Questo script testa tutte le funzionalità della funzione fit_hsmm
# con dati generati random seguendo la logica delle funzioni
# ============================================================================

# Carica il pacchetto R (assumendo che combined_scripts.R sia stato sourcato)
# source("combined_scripts.R")

library(R6)

files <- list.files(path = "R", 
                    pattern = "\\.[Rr]$", 
                    full.names = TRUE,
                    recursive = TRUE)
for (f in files) {
  cat("Importando script:", f, "\n")
  source(f)
}

# Carica librerie necessarie
library(R6)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)
library(tidyr)
library(nnet)
library(MASS)

set.seed(42)

# ============================================================================
# FUNZIONI DI SUPPORTO PER GENERAZIONE DATI
# ============================================================================

#' Genera dati simulati da un HSMM
#' @param n_obs Numero di osservazioni
#' @param n_states Numero di stati
#' @param data_type Tipo di dati ("gaussian", "torus", "poisson", "exponential", "gamma", "beta", "student")
#' @param n_dim Dimensione dei dati (per modelli multivariati)
#' @return Lista con data, true_states, e parametri usati
generate_hsmm_data <- function(n_obs = 200, 
                                n_states = 3, 
                                data_type = "gaussian",
                                n_dim = 2) {
  
  # Genera sequenza di stati con durate variabili
  states <- numeric(n_obs)
  current_state <- sample(1:n_states, 1)
  i <- 1
  
  while(i <= n_obs) {
    # Durata nello stato corrente (esponenziale troncata)
    dwell_time <- min(rpois(1, lambda = 10) + 1, n_obs - i + 1)
    states[i:min(i + dwell_time - 1, n_obs)] <- current_state
    i <- i + dwell_time
    # Transizione a un nuovo stato (non lo stesso)
    current_state <- sample(setdiff(1:n_states, current_state), 1)
  }
  
  # Genera osservazioni basate sugli stati
  if (data_type == "gaussian") {
    # Parametri gaussiani per ogni stato
    mus <- matrix(rnorm(n_states * n_dim, mean = 0, sd = 5), 
                  nrow = n_states, ncol = n_dim)
    data <- matrix(NA, nrow = n_obs, ncol = n_dim)
    for (t in 1:n_obs) {
      data[t, ] <- mvrnorm(1, mu = mus[states[t], ], 
                           Sigma = diag(n_dim))
    }
    
  } else if (data_type == "torus") {
    # Dati circolari (angoli tra 0 e 2*pi)
    mus <- matrix(runif(n_states * 2, 0, 2 * pi), 
                  nrow = n_states, ncol = 2)
    data <- matrix(NA, nrow = n_obs, ncol = 2)
    for (t in 1:n_obs) {
      kappa <- 2  # Concentrazione
      data[t, 1] <- (mus[states[t], 1] + rvonmises(1, mu = 0, kappa = kappa)) %% (2 * pi)
      data[t, 2] <- (mus[states[t], 2] + rvonmises(1, mu = 0, kappa = kappa)) %% (2 * pi)
    }
    
  } else if (data_type == "poisson") {
    # Dati di conteggio
    lambdas <- runif(n_states, min = 2, max = 20)
    data <- matrix(NA, nrow = n_obs, ncol = 1)
    for (t in 1:n_obs) {
      data[t, 1] <- rpois(1, lambda = lambdas[states[t]])
    }
    
  } else if (data_type == "exponential") {
    # Dati positivi con distribuzione esponenziale
    rates <- runif(n_states, min = 0.5, max = 3)
    data <- matrix(NA, nrow = n_obs, ncol = 1)
    for (t in 1:n_obs) {
      data[t, 1] <- rexp(1, rate = rates[states[t]])
    }
    
  } else if (data_type == "gamma") {
    # Dati positivi con distribuzione gamma
    shapes <- runif(n_states, min = 2, max = 5)
    scales <- runif(n_states, min = 0.5, max = 2)
    data <- matrix(NA, nrow = n_obs, ncol = 1)
    for (t in 1:n_obs) {
      data[t, 1] <- rgamma(1, shape = shapes[states[t]], 
                           scale = scales[states[t]])
    }
    
  } else if (data_type == "beta") {
    # Dati in [0, 1]
    alphas <- runif(n_states, min = 2, max = 5)
    betas <- runif(n_states, min = 2, max = 5)
    data <- matrix(NA, nrow = n_obs, ncol = 1)
    for (t in 1:n_obs) {
      data[t, 1] <- rbeta(1, shape1 = alphas[states[t]], 
                          shape2 = betas[states[t]])
    }
    
  } else if (data_type == "student") {
    # Dati con code pesanti (t di Student)
    mus <- rnorm(n_states, mean = 0, sd = 5)
    sigmas <- runif(n_states, min = 0.5, max = 2)
    df <- 5
    data <- matrix(NA, nrow = n_obs, ncol = 1)
    for (t in 1:n_obs) {
      data[t, 1] <- mus[states[t]] + sigmas[states[t]] * rt(1, df = df)
    }
  }
  
  return(list(data = data, true_states = states))
}

# Funzione per von Mises (per dati circolari)
rvonmises <- function(n, mu = 0, kappa = 1) {
  if(kappa < 1e-6) return(runif(n, 0, 2*pi))
  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a)) / (2 * kappa)
  r <- (1 + b^2) / (2 * b)
  
  theta <- numeric(n)
  for(i in 1:n) {
    repeat {
      u1 <- runif(1)
      z <- cos(pi * u1)
      f <- (1 + r * z) / (r + z)
      c <- kappa * (r - f)
      u2 <- runif(1)
      if(c * (2 - c) - u2 > 0) break
      if(log(c / u2) + 1 - c >= 0) break
    }
    u3 <- runif(1)
    theta[i] <- (mu + sign(u3 - 0.5) * acos(f)) %% (2 * pi)
  }
  return(theta)
}

# ============================================================================
# TEST 1: HMM Gaussiano base (semi = FALSE)
# ============================================================================
cat("\n=== TEST 1: HMM Gaussiano base (semi = FALSE) ===\n")

# Genera dati
data1 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                            data_type = "gaussian", n_dim = 2)

# Fit del modello
tryCatch({
  fit1 <- fit_hsmm(
    data = data1$data,
    n_states = 2,
    semi = FALSE,
    max_dwell = NULL,
    model_type = "gaussian",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 1 COMPLETATO CON SUCCESSO\n")
  cat("Log-likelihood finale:", tail(fit1$input_params$llk, 1), "\n")
}, error = function(e) {
  cat("✗ TEST 1 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 2: HSMM Gaussiano base (semi = TRUE, geometric)
# ============================================================================
cat("\n=== TEST 2: HSMM Gaussiano base (semi = TRUE, geometric) ===\n")

# Genera dati
data2 <- generate_hsmm_data(n_obs = 150, n_states = 3, 
                            data_type = "gaussian", n_dim = 2)

# Fit del modello
tryCatch({
  fit2 <- fit_hsmm(
    data = data2$data,
    n_states = 3,
    semi = TRUE,
    max_dwell = 8,
    model_type = "gaussian",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 2 COMPLETATO CON SUCCESSO\n")
  cat("Beta intercepts:", fit2$duration.model$beta.intercepts, "\n")
  cat("Beta time effects:", fit2$duration.model$beta.time_effects, "\n")
}, error = function(e) {
  cat("✗ TEST 2 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 3: HSMM Gaussiano con durate Gompertz (covariates_q)
# ============================================================================
cat("\n=== TEST 3: HSMM Gaussiano con durate Gompertz (covariates_q) ===\n")

# Genera dati con covariate per durate
data3 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                            data_type = "gaussian", n_dim = 2)
covariates_q3 <- matrix(rnorm(150 * 2), nrow = 150, ncol = 2)

# Fit del modello
tryCatch({
  fit3 <- fit_hsmm(
    data = data3$data,
    n_states = 2,
    covariates_q = covariates_q3,
    semi = TRUE,
    max_dwell = 6,
    model_type = "gaussian",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 3 COMPLETATO CON SUCCESSO\n")
  cat("Beta intercepts:", fit3$duration.model$beta.intercepts, "\n")
  cat("Beta covariate:\n")
  print(fit3$duration.model$beta.covariate)
}, error = function(e) {
  cat("✗ TEST 3 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 4: HSMM Gaussiano con transizioni time-varying (covariates_omega)
# ============================================================================
cat("\n=== TEST 4: HSMM Gaussiano con transizioni time-varying (covariates_omega) ===\n")

# Genera dati con covariate per transizioni
data4 <- generate_hsmm_data(n_obs = 150, n_states = 3, 
                            data_type = "gaussian", n_dim = 2)
covariates_omega4 <- rnorm(150)

# Fit del modello
tryCatch({
  fit4 <- fit_hsmm(
    data = data4$data,
    n_states = 3,
    covariates_omega = covariates_omega4,
    semi = TRUE,
    max_dwell = 7,
    model_type = "gaussian",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 4 COMPLETATO CON SUCCESSO\n")
  cat("Dimensioni omega array:", dim(fit4$transition.model$omega), "\n")
}, error = function(e) {
  cat("✗ TEST 4 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 5: HSMM Gaussiano completo (entrambe le covariate)
# ============================================================================
cat("\n=== TEST 5: HSMM Gaussiano completo (entrambe le covariate) ===\n")

# Genera dati con entrambe le covariate
data5 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                            data_type = "gaussian", n_dim = 2)
covariates_omega5 <- rnorm(150)
covariates_q5 <- matrix(rnorm(150 * 1), nrow = 150, ncol = 1)

# Fit del modello
tryCatch({
  fit5 <- fit_hsmm(
    data = data5$data,
    n_states = 2,
    covariates_omega = covariates_omega5,
    covariates_q = covariates_q5,
    semi = TRUE,
    max_dwell = 6,
    model_type = "gaussian",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 5 COMPLETATO CON SUCCESSO\n")
  cat("Beta covariate:\n")
  print(fit5$duration.model$beta.covariate)
  cat("Dimensioni omega array:", dim(fit5$transition.model$omega), "\n")
}, error = function(e) {
  cat("✗ TEST 5 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 6: HSMM Torus (dati circolari)
# ============================================================================
cat("\n=== TEST 6: HSMM Torus (dati circolari) ===\n")

# Genera dati circolari
data6 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                            data_type = "torus")

# Fit del modello
tryCatch({
  fit6 <- fit_hsmm(
    data = data6$data,
    n_states = 2,
    semi = TRUE,
    max_dwell = 8,
    model_type = "torus",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 6 COMPLETATO CON SUCCESSO\n")
  cat("Tau mu1:", fit6$emission.model$tau.mu1, "\n")
  cat("Tau k1:", fit6$emission.model$tau.k1, "\n")
}, error = function(e) {
  cat("✗ TEST 6 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 7: HSMM Poisson (dati di conteggio)
# ============================================================================
cat("\n=== TEST 7: HSMM Poisson (dati di conteggio) ===\n")

# Genera dati Poisson
data7 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                            data_type = "poisson")

# Fit del modello
tryCatch({
  fit7 <- fit_hsmm(
    data = data7$data,
    n_states = 2,
    semi = TRUE,
    max_dwell = 7,
    model_type = "poisson",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 7 COMPLETATO CON SUCCESSO\n")
  cat("Lambda stimati:", fit7$emission.model$lambda, "\n")
}, error = function(e) {
  cat("✗ TEST 7 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 8: HSMM Exponential (dati positivi)
# ============================================================================
cat("\n=== TEST 8: HSMM Exponential (dati positivi) ===\n")

# Genera dati esponenziali
data8 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                            data_type = "exponential")

# Fit del modello
tryCatch({
  fit8 <- fit_hsmm(
    data = data8$data,
    n_states = 2,
    semi = TRUE,
    max_dwell = 6,
    model_type = "exponential",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 8 COMPLETATO CON SUCCESSO\n")
  cat("Lambda stimati:", fit8$emission.model$lambda, "\n")
}, error = function(e) {
  cat("✗ TEST 8 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 9: HSMM Gamma (dati positivi)
# ============================================================================
cat("\n=== TEST 9: HSMM Gamma (dati positivi) ===\n")

# Genera dati gamma
data9 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                            data_type = "gamma")

# Fit del modello
tryCatch({
  fit9 <- fit_hsmm(
    data = data9$data,
    n_states = 2,
    semi = TRUE,
    max_dwell = 7,
    model_type = "gamma",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 9 COMPLETATO CON SUCCESSO\n")
  cat("Shape stimati:", fit9$emission.model$shape, "\n")
  cat("Scale stimati:", fit9$emission.model$scale, "\n")
}, error = function(e) {
  cat("✗ TEST 9 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 10: HSMM Beta (dati in [0,1])
# ============================================================================
cat("\n=== TEST 10: HSMM Beta (dati in [0,1]) ===\n")

# Genera dati beta
data10 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                             data_type = "beta")

# Fit del modello
tryCatch({
  fit10 <- fit_hsmm(
    data = data10$data,
    n_states = 2,
    semi = TRUE,
    max_dwell = 6,
    model_type = "beta",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 10 COMPLETATO CON SUCCESSO\n")
  cat("Alpha stimati:", fit10$emission.model$alpha, "\n")
  cat("Beta stimati:", fit10$emission.model$beta, "\n")
}, error = function(e) {
  cat("✗ TEST 10 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST 11: HSMM Student-t (dati con code pesanti)
# ============================================================================
cat("\n=== TEST 11: HSMM Student-t (dati con code pesanti) ===\n")

# Genera dati Student-t
data11 <- generate_hsmm_data(n_obs = 150, n_states = 2, 
                             data_type = "student")

# Fit del modello
tryCatch({
  fit11 <- fit_hsmm(
    data = data11$data,
    n_states = 2,
    semi = TRUE,
    max_dwell = 7,
    model_type = "student",
    max_iter = 20,
    tol = 1e-4,
    verbose = TRUE,
    seed = 123
  )
  cat("✓ TEST 11 COMPLETATO CON SUCCESSO\n")
  cat("Mu stimati:", fit11$emission.model$mu, "\n")
  cat("Sigma stimati:", fit11$emission.model$sigma, "\n")
  cat("Df:", fit11$emission.model$df, "\n")
}, error = function(e) {
  cat("✗ TEST 11 FALLITO:", e$message, "\n")
})

# ============================================================================
# TEST AGGIUNTIVI: Verifica struttura output
# ============================================================================
cat("\n=== TEST STRUTTURA OUTPUT ===\n")

# Usa fit2 come esempio
if(exists("fit2")) {
  cat("\nStruttura dell'output di fit_hsmm:\n")
  cat("Componenti principali:", names(fit2), "\n\n")
  
  cat("input_params contiene:", names(fit2$input_params), "\n\n")
  
  cat("emission.model è di classe:", class(fit2$emission.model), "\n")
  cat("emission.model contiene:", names(fit2$emission.model), "\n\n")
  
  cat("transition.model è di classe:", class(fit2$transition.model), "\n")
  cat("transition.model contiene:", names(fit2$transition.model), "\n\n")
  
  cat("duration.model è di classe:", class(fit2$duration.model), "\n")
  cat("duration.model contiene:", names(fit2$duration.model), "\n\n")
  
  cat("Dimensioni della matrice gamma:", dim(fit2$transition.model$gamma), "\n")
  cat("Dimensioni di omega:", 
      if(is.matrix(fit2$transition.model$omega)) {
        dim(fit2$transition.model$omega)
      } else {
        "omega è un array"
      }, "\n")
  cat("Dimensioni di post.pi:", dim(fit2$transition.model$post.pi), "\n")
  cat("Lunghezza di Pi:", length(fit2$transition.model$Pi), "\n")
}

cat("\n=== TUTTI I TEST COMPLETATI ===\n")
cat("Script di test per fit_hsmm terminato.\n")
cat("Verifica i messaggi sopra per vedere quali test sono stati completati con successo.\n")
