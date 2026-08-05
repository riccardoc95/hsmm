risultati <- lapply(1:100, function(b) {

  sim_b <- simulate_hsmm(
    n = 500,
    n_states = 3,
    model_type = "gaussian",
    max_dwell = 12,
    seed = b
  )

  fit_b <- fit_hsmm(
    sim_b$data,
    n_states = 3,
    model_type = "gaussian",
    max_dwell = 12,
    max_iter = 200,
    tol = 1e-6,
    verbose = FALSE,
    seed = 1000 + b
  )

  tab <- table(
    factor(sim_b$state, levels = 1:3),
    factor(fit_b$decoded_state, levels = 1:3)
  )

  perm <- clue::solve_LSAP(tab, maximum = TRUE)

  mappa <- integer(3)
  for (k in 1:3) {
    mappa[perm[k]] <- k
  }

  decoded <- mappa[fit_b$decoded_state]

  c(
    converged = fit_b$converged,
    accuracy = mean(decoded == sim_b$state),
    loglik = fit_b$final_loglik
  )
})

risultati <- do.call(rbind, risultati)

colMeans(risultati)
summary(risultati[, "accuracy"])
mean(risultati[, "converged"])
