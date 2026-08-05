expectation.step <- function(params, emission.model, transition.model) {
  fit <- emission.model$density
  fit[!is.finite(fit) | fit <= 0] <- .Machine$double.xmin

  bw <- backward.forward(
    fit = fit,
    Gamma = transition.model$gamma,
    Pi = transition.model$Pi,
    K = params$n_states,
    M = params$max_dwell
  )
  bw$post.pi <- .collapse_posterior(
    bw$post.pi.aug, params$n_states, params$max_dwell
  )
  bw
}

em.step <- function(params, emission.model, transition.model, duration.model) {
  bw <- expectation.step(params, emission.model, transition.model)

  duration.model$compute_beta(params, bw$post.bi.pi.aug)
  duration.model$compute_p.array(params)

  transition.model$compute_omega(params, bw$post.bi.pi.aug)
  transition.model$compute_gamma(params, duration.model$p.array)
  transition.model$update(params, bw$post.pi.aug)

  emission.model$update(params, transition.model$Pi, bw$post.pi)
  bw
}

em.algorithm <- function(params, emission.model, transition.model, duration.model) {
  llk <- numeric(params$max_iter)
  converged <- FALSE
  used <- params$max_iter

  for (step in seq_len(params$max_iter)) {
    bw <- em.step(params, emission.model, transition.model, duration.model)
    llk[step] <- bw$l

    dif <- if (step == 1) NA_real_ else llk[step] - llk[step - 1]
    if (params$verbose) {
      message <- paste("iteration", step, "; loglik =", signif(llk[step], 8))
      if (step > 1) message <- paste(message, "; diff =", signif(dif, 5))
      cat(message, "\n")
    }

    if (step > 1 && is.finite(dif) && abs(dif) < params$tol) {
      converged <- TRUE
      used <- step
      break
    }
  }

  final <- expectation.step(params, emission.model, transition.model)
  transition.model$post.pi <- final$post.pi

  list(
    loglik = llk[seq_len(used)],
    final_loglik = final$l,
    converged = converged,
    iterations = used,
    posterior = final$post.pi,
    posterior_augmented = final$post.pi.aug,
    pairwise_posterior = final$post.bi.pi.aug
  )
}
