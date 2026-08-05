sim_covariates <- function(x, n) {
  if (is.null(x)) return(matrix(numeric(0), n, 0))
  x <- .as_numeric_matrix(x)
  if (nrow(x) == n - 1) x <- rbind(x, x[nrow(x), , drop = FALSE])
  if (nrow(x) != n) stop("covariates must have n or n - 1 rows")
  x
}

sim_default_emission <- function(model_type, n_states, n_dim, ar) {
  centers <- if (n_states == 1) 0 else seq(-3, 3, length.out = n_states)
  mu <- matrix(rep(centers, n_dim), n_states, n_dim)
  for (j in seq_len(n_dim)) mu[, j] <- mu[, j] + (j - 1) / 2

  if (model_type == "gaussian") {
    return(list(mu = mu, sigma = rep(list(diag(n_dim) * 0.7), n_states)))
  }

  if (model_type == "poisson") {
    lambda <- matrix(rep(seq(2, 12, length.out = n_states), n_dim),
                     n_states, n_dim)
    return(list(lambda = lambda))
  }

  if (model_type == "exponential") {
    rate <- matrix(rep(seq(1.8, 0.35, length.out = n_states), n_dim),
                   n_states, n_dim)
    return(list(rate = rate))
  }

  if (model_type == "gamma") {
    shape <- matrix(4, n_states, n_dim)
    means <- matrix(rep(seq(1, 8, length.out = n_states), n_dim),
                    n_states, n_dim)
    return(list(shape = shape, rate = shape / means))
  }

  if (model_type == "beta") {
    means <- matrix(rep(seq(0.2, 0.8, length.out = n_states), n_dim),
                    n_states, n_dim)
    return(list(shape1 = means * 20, shape2 = (1 - means) * 20))
  }

  if (model_type == "student") {
    return(list(mu = mu,
                scale = matrix(0.8, n_states, n_dim),
                df = rep(5, n_states)))
  }

  if (model_type == "torus") {
    base_mu <- seq(0, 2 * pi, length.out = n_states + 1)[seq_len(n_states)]
    for (j in seq_len(n_dim)) {
      mu[, j] <- (base_mu + (j - 1) * pi / 5) %% (2 * pi)
    }
    return(list(mu = mu, kappa = matrix(8, n_states, n_dim)))
  }

  coef <- vector("list", n_states)
  sigma <- rep(list(diag(n_dim) * 0.5), n_states)
  for (k in seq_len(n_states)) {
    B <- matrix(0, 1 + n_dim * ar, n_dim)
    B[1, ] <- centers[k] * 0.35
    for (j in seq_len(n_dim)) B[1 + j, j] <- 0.65
    coef[[k]] <- B
  }
  list(coef = coef, sigma = sigma)
}

#' Simulate data from an HSMM
#'
#' @export
simulate_hsmm <- function(
    n = 300,
    n_states = 3,
    model_type = "gaussian",
    n_dim = 1,
    semi = TRUE,
    max_dwell = 20,
    covariates_omega = NULL,
    covariates_q = NULL,
    Pi = NULL,
    omega = NULL,
    duration_intercept = NULL,
    duration_age_effect = NULL,
    duration_coef = NULL,
    transition_coef = NULL,
    emission = NULL,
    ar = 1,
    seed = NULL,
    family = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (!is.null(family)) model_type <- family

  model_type <- tolower(model_type)

  allowed <- c("gaussian", "poisson", "exponential", "gamma",
               "beta", "student", "torus")#, "vargaussian")
  if (!model_type %in% allowed) stop("unknown model_type")

  n <- as.integer(n)
  n_states <- as.integer(n_states)
  n_dim <- as.integer(n_dim)
  ar <- as.integer(ar)
  max_dwell <- as.integer(max_dwell)
  if (!semi) max_dwell <- 1L

  # if (model_type == "vargaussian" && ar < 1) stop("ar must be positive")

  covariates_q <- sim_covariates(covariates_q, n)
  covariates_omega <- sim_covariates(covariates_omega, n)

  if (is.null(Pi)) Pi <- rep(1 / n_states, n_states)
  Pi <- .normalize_probability(Pi)

  if (is.null(omega)) {
    if (n_states == 1) {
      omega <- matrix(1, 1, 1)
    } else {
      omega <- matrix(1 / (n_states - 1), n_states, n_states)
      diag(omega) <- 0
    }
  }
  omega <- .normalize_omega(omega)

  if (is.null(duration_intercept)) duration_intercept <- rep(-2.1, n_states)
  if (is.null(duration_age_effect)) {
    duration_age_effect <- if (semi) rep(0.12, n_states) else rep(0, n_states)
  }
  if (is.null(duration_coef)) {
    duration_coef <- matrix(0, n_states, ncol(covariates_q))
  } else {
    duration_coef <- as.matrix(duration_coef)
  }

  default_emission <- sim_default_emission(model_type, n_states, n_dim, ar)
  if (!is.null(emission)) {
    for (name in names(emission)) default_emission[[name]] <- emission[[name]]
  }
  emission <- default_emission

  state <- integer(n)
  dwell_age <- integer(n)
  state[1] <- sample.int(n_states, 1, prob = Pi)
  dwell_age[1] <- 1L

  for (t in 2:n) {
    k <- state[t - 1]
    age <- min(dwell_age[t - 1], max_dwell)

    eta <- duration_intercept[k] + duration_age_effect[k] * (age + 0.5)
    if (ncol(covariates_q) > 0) {
      eta <- eta + sum(covariates_q[t - 1, ] * duration_coef[k, ])
    }
    eta <- max(min(eta, 20), -20)
    p_exit <- 1 - exp(-exp(eta))

    if (n_states == 1) {
      state_probability <- 1
    } else {
      destinations <- setdiff(seq_len(n_states), k)
      destination_probability <- omega[k, ]

      if (!is.null(transition_coef) && length(destinations) > 1) {
        x <- c(1, covariates_omega[t - 1, ])
        score <- c(drop(x %*% transition_coef[[k]]), 0)
        score <- exp(score - max(score))
        destination_probability[] <- 0
        destination_probability[destinations] <- score / sum(score)
      }

      state_probability <- p_exit * destination_probability
      state_probability[k] <- 1 - p_exit
    }

    state[t] <- sample.int(n_states, 1, prob = state_probability)
    if (state[t] == k) {
      dwell_age[t] <- min(dwell_age[t - 1] + 1L, max_dwell)
    } else {
      dwell_age[t] <- 1L
    }
  }

  data <- matrix(NA_real_, n, n_dim)
  for (t in seq_len(n)) {
    k <- state[t]

    if (model_type == "gaussian") {
      data[t, ] <- .rmvnorm_basic(1, emission$mu[k, ], emission$sigma[[k]])
    } else if (model_type == "poisson") {
      data[t, ] <- stats::rpois(n_dim, emission$lambda[k, ])
    } else if (model_type == "exponential") {
      data[t, ] <- stats::rexp(n_dim, emission$rate[k, ])
    } else if (model_type == "gamma") {
      data[t, ] <- stats::rgamma(n_dim, emission$shape[k, ], emission$rate[k, ])
    } else if (model_type == "beta") {
      data[t, ] <- stats::rbeta(n_dim, emission$shape1[k, ], emission$shape2[k, ])
    } else if (model_type == "student") {
      data[t, ] <- emission$mu[k, ] +
        emission$scale[k, ] * stats::rt(n_dim, emission$df[k])
    } else if (model_type == "torus") {
      for (j in seq_len(n_dim)) {
        data[t, j] <- .rvonmises_basic(1, emission$mu[k, j], emission$kappa[k, j])
      }
    } else {
      if (t <= ar) {
        mean_t <- emission$coef[[k]][1, ]
      } else {
        x <- 1
        for (lag in seq_len(ar)) x <- c(x, data[t - lag, ])
        mean_t <- drop(x %*% emission$coef[[k]])
      }
      data[t, ] <- .rmvnorm_basic(1, mean_t, emission$sigma[[k]])
    }
  }
  colnames(data) <- paste0("y", seq_len(n_dim))

  out <- list(
    data = data,
    state = state,
    dwell_age = dwell_age,
    covariates_omega = covariates_omega,
    covariates_q = covariates_q,
    parameters = list(
      n_states = n_states,
      model_type = model_type,
      n_dim = n_dim,
      semi = isTRUE(semi),
      max_dwell = max_dwell,
      Pi = Pi,
      omega = omega,
      duration_intercept = duration_intercept,
      duration_age_effect = duration_age_effect,
      duration_coef = duration_coef,
      transition_coef = transition_coef,
      emission = emission,
      ar = ar
    )
  )
  class(out) <- "hsmm_simulation"
  out
}

sim_emission_from_fit <- function(object) {
  model <- object$emission.model
  switch(
    object$input_params$model_type,
    gaussian = list(mu = model$mu, sigma = model$sigma),
    poisson = list(lambda = model$lambda),
    exponential = list(rate = model$rate),
    gamma = list(shape = model$shape, rate = model$rate),
    beta = list(shape1 = model$shape1, shape2 = model$shape2),
    student = list(mu = model$mu, scale = model$scale, df = model$df),
    torus = list(mu = model$mu, kappa = model$kappa),
    # vargaussian = list(coef = model$coef, sigma = model$sigma)
  )
}

.sim_reuse_covariates <- function(new, old, n) {
  if (!is.null(new) || is.null(old) || ncol(old) == 0) return(new)
  old[rep(seq_len(nrow(old)), length.out = n), , drop = FALSE]
}

#' Simulate from a fitted HSMM
#'
#' @export
simulate.hsmm_fit <- function(object, nsim = 1, seed = NULL,
                              n = object$input_params$n_obs,
                              covariates_omega = NULL,
                              covariates_q = NULL, ...) {
  nsim <- as.integer(nsim)
  if (nsim < 1) stop("nsim must be positive")

  covariates_omega <- .sim_reuse_covariates(
    covariates_omega, object$input_params$covariates_omega, n
  )
  covariates_q <- .sim_reuse_covariates(
    covariates_q, object$input_params$covariates_q, n
  )

  simulations <- vector("list", nsim)
  for (i in seq_len(nsim)) {
    new_seed <- if (is.null(seed)) NULL else seed + i - 1
    simulations[[i]] <- simulate_hsmm(
      n = n,
      n_states = object$input_params$n_states,
      model_type = object$input_params$model_type,
      n_dim = ncol(object$input_params$data),
      semi = object$input_params$semi,
      max_dwell = object$input_params$max_dwell,
      covariates_omega = covariates_omega,
      covariates_q = covariates_q,
      Pi = object$transition.model$Pi,
      omega = object$transition.model$omega,
      duration_intercept = object$duration.model$beta.intercepts,
      duration_age_effect = object$duration.model$beta.time_effects,
      duration_coef = object$duration.model$beta.covariate,
      transition_coef = object$transition.model$coefficients,
      emission = sim_emission_from_fit(object),
      ar = object$input_params$ar,
      seed = new_seed
    )
  }

  if (nsim == 1) simulations[[1]] else simulations
}

#' @export
print.hsmm_simulation <- function(x, ...) {
  cat("Simulated hidden semi-Markov series\n")
  cat("  observations:", nrow(x$data), "\n")
  cat("  states:", x$parameters$n_states, "\n")
  cat("  emission:", x$parameters$model_type, "\n")
  invisible(x)
}

#' @export
plot.hsmm_simulation <- function(x, ...) {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))
  graphics::par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
  graphics::plot(x$data[, 1], type = "l", xlab = "time",
                 ylab = "observation", main = "Simulated observations", ...)
  graphics::plot(x$state, type = "s", xlab = "time",
                 ylab = "state", main = "Simulated latent state")
  invisible(x)
}
