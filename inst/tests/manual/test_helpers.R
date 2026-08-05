# Small helper functions used by the manual tests.

make_permutations <- function(x) {
  if (length(x) == 1) {
    return(matrix(x, nrow = 1))
  }

  out <- NULL

  for (i in seq_along(x)) {
    smaller <- make_permutations(x[-i])
    first_column <- rep(x[i], nrow(smaller))
    out <- rbind(out, cbind(first_column, smaller))
  }

  out
}

align_states <- function(real_state, fitted_state, n_states) {
  permutations <- make_permutations(seq_len(n_states))

  best_accuracy <- -1
  best_state <- fitted_state
  best_permutation <- seq_len(n_states)

  for (i in seq_len(nrow(permutations))) {
    new_state <- fitted_state

    for (k in seq_len(n_states)) {
      new_state[fitted_state == k] <- permutations[i, k]
    }

    accuracy <- mean(real_state == new_state)

    if (accuracy > best_accuracy) {
      best_accuracy <- accuracy
      best_state <- new_state
      best_permutation <- permutations[i, ]
    }
  }

  list(
    accuracy = best_accuracy,
    state = best_state,
    permutation = best_permutation
  )
}

check_fit <- function(sim, fit, minimum_accuracy = 0.75,
                      test_name = "HSMM test", make_plot = TRUE) {
  n_states <- sim$parameters$n_states
  aligned <- align_states(sim$state, fit$decoded_state, n_states)

  posterior_error <- max(abs(rowSums(fit$posterior) - 1))
  finite_loglik <- is.finite(fit$final_loglik)
  loglik_gain <- fit$final_loglik - fit$loglik[1]

  if (length(fit$loglik) > 1) {
    smallest_loglik_step <- min(diff(fit$loglik))
  } else {
    smallest_loglik_step <- NA_real_
  }

  passed <- finite_loglik &&
    posterior_error < 1e-6 &&
    aligned$accuracy >= minimum_accuracy

  cat("\n========================================\n")
  cat(test_name, "\n")
  cat("States:", n_states, "\n")
  cat("Observed variables:", ncol(sim$data), "\n")
  cat("State accuracy:", round(aligned$accuracy, 3), "\n")
  cat("Posterior error:", signif(posterior_error, 3), "\n")
  cat("Final log-likelihood:", round(fit$final_loglik, 3), "\n")
  cat("Log-likelihood gain:", round(loglik_gain, 3), "\n")
  cat("Smallest EM step:", round(smallest_loglik_step, 6), "\n")
  cat("Converged:", fit$converged, "\n")
  cat("Iterations:", fit$iterations, "\n")
  cat("Result:", if (passed) "OK" else "NOT OK", "\n")

  print(table(
    simulated = sim$state,
    fitted = aligned$state
  ))

  if (make_plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

    plot(
      sim$data[, 1],
      type = "l",
      main = "First observed variable",
      xlab = "time",
      ylab = "observation"
    )

    plot(
      sim$state,
      type = "s",
      ylim = c(0.8, n_states + 0.2),
      main = "Simulated latent state",
      xlab = "time",
      ylab = "state"
    )

    plot(
      aligned$state,
      type = "s",
      ylim = c(0.8, n_states + 0.2),
      main = "Decoded latent state",
      xlab = "time",
      ylab = "state"
    )
  }

  result <- list(
    passed = passed,
    accuracy = aligned$accuracy,
    aligned_state = aligned$state,
    permutation = aligned$permutation,
    posterior_error = posterior_error,
    loglik_gain = loglik_gain,
    smallest_loglik_step = smallest_loglik_step,
    simulation = sim,
    fit = fit
  )

  if (!passed) {
    stop(test_name, " did not pass")
  }

  invisible(result)
}
