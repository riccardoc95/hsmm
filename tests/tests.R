library(R6)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)
library(tidyr)
library(nnet)
library(MASS)
library(circular)
library(jsonlite)
library(hsmm)
#
# devtools::document()
# roxygen2::roxygenise()
# devtools::build_manual()
# R CMD Rd2pdf .
# install.packages("tinytex")
# tinytex::install_tinytex()
# devtools::check(manual = TRUE)
# pkgdown::build_site()

# usethis::create_package("~/Develops/rstudio/hsmm/hsmm")
# devtools::document()
# devtools::check()
# devtools::build()
# install.packages("hsmm_0.1.0.tar.gz", repos = NULL, type = "source")

# devtools::install()

# setwd("~/Develops/rstudio/hsmm/hsmm")
#
# set.seed(42)
#
# files <- list.files(path = "R",
#                     pattern = "\\.[Rr]$",
#                     full.names = TRUE,
#                     recursive = TRUE)
# for (f in files) {
#   cat("Importando script:", f, "\n")
#   source(f)
# }


generate_data <- function(n_obs = 200,
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

  return(data)
}



test_fit_hsmm_errors <- function(model_type="torus") {
  set.seed(42)
  n_obs <- 50
  n_states_list <- c(2, 3, 8)  # variazione di n_states

  covariates_omega_list <- list(NULL,
                                matrix(runif(n_obs*2), ncol=2),
                                matrix(runif(n_obs*4), ncol=4))
  covariates_q_list <- list(NULL,
                            matrix(runif(n_obs*2), ncol=2),
                            matrix(runif(n_obs*4), ncol=4))
  semi_list <- c(TRUE, FALSE)
  max_dwell_list <- list(2, 5, 8)
  max_iter_list <- c(200)
  tol_list <- c(1e-5)
  verbose_list <- c(TRUE)
  seeds <- c(1234)

  error_log <- list()
  success_log <- list()
  comb_num <- 0

  for (n_states in n_states_list) {
    # Rigenera dati ad angoli per il numero di stati corrente
    data <- generate_data(n_obs=n_obs,  n_states = n_states,
                  data_type = model_type, n_dim = 2)
    #n_features <- 2  # o numero di dimensioni volute
    #data <- matrix(rnorm(n_obs * n_features), ncol = n_features)
    for (cov_omega in covariates_omega_list) {
      for (cov_q in covariates_q_list) {
        for (semi_val in semi_list) {
          for (max_dw in max_dwell_list) {
            for (max_it in max_iter_list) {
              for (tol_val in tol_list) {
                for (verb in verbose_list) {
                  for (sd in seeds) {

                    comb_num <- comb_num + 1
                    cat(sprintf("Testing combo #%d with n_states=%d...\n", comb_num, n_states))

                    try_res <- tryCatch({
                      fit_res <- fit_hsmm(data = data,
                                          n_states = n_states,
                                          covariates_omega = cov_omega,
                                          covariates_q = cov_q,
                                          semi = semi_val,
                                          max_dwell = max_dw,
                                          model_type = model_type,
                                          max_iter = max_it,
                                          tol = tol_val,
                                          verbose = verb,
                                          seed = sd)
                      list(success=TRUE, result=fit_res)
                    }, error = function(e) {
                      calls <- sys.calls()
                      calls_text <- capture.output(print(calls))
                      list(success=FALSE,
                           message = e$message,
                           traceback = calls_text,
                           combination = list(
                             n_states = n_states,
                             covariates_omega = cov_omega,
                             covariates_q = cov_q,
                             semi = semi_val,
                             max_dwell = max_dw,
                             max_iter = max_it,
                             tol = tol_val,
                             verbose = verb,
                             seed = sd
                           ))
                    })

                    if (!try_res$success) {
                      error_log[[length(error_log) + 1]] <- list(
                        combo_number = comb_num,
                        error_message = try_res$message,
                        traceback = try_res$traceback,
                        params = list(
                          data = data,
                          n_states = n_states,
                          covariates_omega = cov_omega,
                          covariates_q = cov_q,
                          semi = semi_val,
                          max_dwell = max_dw,
                          max_iter = max_it,
                          model_type = model_type,
                          tol = tol_val,
                          verbose = verb,
                          seed = sd
                        )
                      )
                      cat(sprintf("Error in combo #%d.\n", comb_num))
                    } else {
                      success_log[[length(success_log) + 1]] <- list(
                        combo_number = comb_num,
                        params = list(
                          data = data,
                          n_states = n_states,
                          covariates_omega = if(is.null(cov_omega)) "NULL" else "matrix",
                          covariates_q = if(is.null(cov_q)) "NULL" else "matrix",
                          semi = semi_val,
                          max_dwell = max_dw,
                          max_iter = max_it,
                          model_type = model_type,
                          tol = tol_val,
                          verbose = verb,
                          seed = sd
                        )#,
                        #result = try_res$result  # opzionale, puoi scegliere cosa salvare
                      )
                      cat(sprintf("Combo #%d success.\n", comb_num))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  cat("\n--- TEST LOG SUMMARY ---\n")
  cat(sprintf("Successes: %d\n", length(success_log)))
  cat(sprintf("Errors: %d\n", length(error_log)))

  # Salvare entrambi i log in un'unica lista da scrivere in JSON
  test_log <- list(
    successes = success_log,
    errors = error_log
  )
  return(test_log)
}

test_log <- test_fit_hsmm_errors(model_type="gaussian")


# write_json(test_log, path = "test_log.json", pretty = TRUE, auto_unbox = TRUE)
