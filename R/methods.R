#' @export
print.hsmm_fit <- function(x, ...) {
  cat("Hidden semi-Markov model\n")
  cat("  observations:", x$input_params$n_obs, "\n")
  cat("  states:", x$input_params$n_states, "\n")
  cat("  emission:", x$input_params$model_type, "\n")
  cat("  max dwell:", x$input_params$max_dwell, "\n")
  cat("  iterations:", x$iterations, "\n")
  cat("  final log-likelihood:", signif(x$final_loglik, 8), "\n")
  invisible(x)
}

#' @export
summary.hsmm_fit <- function(object, ...) {
  answer <- list(
    call = object$call,
    n_obs = object$input_params$n_obs,
    n_states = object$input_params$n_states,
    emission = object$input_params$model_type,
    max_dwell = object$input_params$max_dwell,
    iterations = object$iterations,
    converged = object$converged,
    loglik = object$final_loglik,
    initial_probability = object$transition.model$Pi,
    transition_probability = object$transition.model$omega,
    duration_intercept = object$duration.model$beta.intercepts,
    duration_age_effect = object$duration.model$beta.time_effects
  )
  class(answer) <- "summary.hsmm_fit"
  answer
}

#' @export
print.summary.hsmm_fit <- function(x, ...) {
  cat("Summary of hsmm_fit\n")
  cat("  observations:", x$n_obs, "\n")
  cat("  states:", x$n_states, "\n")
  cat("  emission:", x$emission, "\n")
  cat("  converged:", x$converged, "\n")
  cat("  log-likelihood:", signif(x$loglik, 8), "\n\n")
  cat("Initial probabilities\n")
  print(x$initial_probability)
  cat("\nMean conditional transition matrix\n")
  print(x$transition_probability)
  invisible(x)
}

#' @export
predict.hsmm_fit <- function(object, type = c("state", "posterior"), ...) {
  type <- match.arg(type)
  if (type == "state") object$decoded_state else object$posterior
}

#' @export
plot.hsmm_fit <- function(x, ...) {
  y <- x$input_params$data[, 1]
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))
  graphics::par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
  graphics::plot(y, type = "l", xlab = "time", ylab = "observation",
                 main = "First observed variable", ...)
  graphics::plot(x$decoded_state, type = "s", xlab = "time", ylab = "state",
                 main = "Decoded latent state")
  invisible(x)
}
