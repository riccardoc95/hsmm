library(hsmm)

n <- 400
n_states <- 2
max_dwell <- 10
max_iter <- 40
minimum_accuracy <- 0.75

families <- c(
  "gaussian",
  "poisson",
  "exponential",
  "gamma",
  "beta",
  "student",
  "torus"
)

make_emission <- function(family) {
  if (family == "gaussian") {
    return(list(
      mu = matrix(c(-4, 4), nrow = 2),
      sigma = list(matrix(0.5, 1, 1), matrix(0.5, 1, 1))
    ))
  }

  if (family == "poisson") {
    return(list(
      lambda = matrix(c(2, 15), nrow = 2)
    ))
  }

  if (family == "exponential") {
    return(list(
      rate = matrix(c(2.5, 0.25), nrow = 2)
    ))
  }

  if (family == "gamma") {
    return(list(
      shape = matrix(c(6, 6), nrow = 2),
      rate = matrix(c(7.5, 0.75), nrow = 2)
    ))
  }

  if (family == "beta") {
    return(list(
      shape1 = matrix(c(4.5, 25.5), nrow = 2),
      shape2 = matrix(c(25.5, 4.5), nrow = 2)
    ))
  }

  if (family == "student") {
    return(list(
      mu = matrix(c(-4, 4), nrow = 2),
      scale = matrix(c(0.6, 0.6), nrow = 2),
      df = c(7, 7)
    ))
  }

  if (family == "torus") {
    return(list(
      mu = matrix(c(pi / 2, 3 * pi / 2), nrow = 2),
      kappa = matrix(c(12, 12), nrow = 2)
    ))
  }
}


state_accuracy <- function(real_state, fitted_state) {
  accuracy_1 <- mean(real_state == fitted_state)
  accuracy_2 <- mean(real_state == (3 - fitted_state))
  max(accuracy_1, accuracy_2)
}

results <- data.frame()
all_runs <- list()

for (i in seq_along(families)) {
  family <- families[i]
  cat("\n----------------------------------------\n")
  cat("Family:", family, "\n")

  one_result <- tryCatch({
    sim <- simulate_hsmm(
      n = n,
      n_states = n_states,
      model_type = family,
      n_dim = 1,
      max_dwell = max_dwell,
      Pi = c(0.5, 0.5),
      duration_intercept = c(-2.4, -2.2),
      duration_age_effect = c(0.18, 0.14),
      emission = make_emission(family),
      seed = 100 + i
    )

    fit <- fit_hsmm(
      data = sim$data,
      n_states = n_states,
      model_type = family,
      max_dwell = max_dwell,
      max_iter = max_iter,
      tol = 1e-4,
      verbose = FALSE,
      seed = 200 + i
    )

    accuracy <- state_accuracy(sim$state, fit$decoded_state)
    posterior_error <- max(abs(rowSums(fit$posterior) - 1))
    loglik_gain <- fit$final_loglik - fit$loglik[1]

    passed <- is.finite(fit$final_loglik) &&
      posterior_error < 1e-6 &&
      accuracy >= minimum_accuracy

    all_runs[[family]] <- list(
      simulation = sim,
      fit = fit,
      accuracy = accuracy
    )

    cat("Accuratezza stati:", round(accuracy, 3), "\n")
    cat("Errore posteriori:", signif(posterior_error, 3), "\n")
    cat("Aumento log-likelihood:", round(loglik_gain, 3), "\n")
    cat("Convergenza:", fit$converged, "\n")
    cat("Risultato:", if (passed) "OK" else "NON OK", "\n")

    data.frame(
      family = family,
      accuracy = round(accuracy, 3),
      posterior_error = posterior_error,
      loglik_gain = round(loglik_gain, 3),
      converged = fit$converged,
      iterations = fit$iterations,
      passed = passed,
      error = ""
    )
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")

    data.frame(
      family = family,
      accuracy = NA_real_,
      posterior_error = NA_real_,
      loglik_gain = NA_real_,
      converged = FALSE,
      iterations = NA_integer_,
      passed = FALSE,
      error = conditionMessage(e)
    )
  })

  results <- rbind(results, one_result)
}

cat("\n\nFINAL RESULTS\n")
print(results, row.names = FALSE)

if (all(results$passed)) {
  cat("\nTutti i casi sono stati superati.\n")
} else {
  stop("Almeno un caso non ha superato il test.")
}
