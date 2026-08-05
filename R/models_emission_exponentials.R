GaussianModel <- R6::R6Class(
  "GaussianModel",
  inherit = EmissionModel,
  public = list(
    mu = NULL,
    sigma = NULL,

    initialize = function(params, Pi, post.pi) {
      K <- params$n_states
      self$mu <- matrix(0, K, ncol(params$data))
      self$sigma <- vector("list", K)
      self$update(params, Pi, post.pi)
      .copy_init_fields(self, params$init$emission)
      self$compute_density(params)
    },

    compute_density = function(params) {
      K <- params$n_states
      self$density <- matrix(0, params$n_obs, K)
      for (k in seq_len(K)) {
        self$density[, k] <- .dmvnorm_basic(
          params$data, self$mu[k, ], self$sigma[[k]]
        )
      }
      self$density <- pmax(self$density, .Machine$double.xmin)
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      for (k in seq_len(params$n_states)) {
        w <- post.pi[, k]
        self$mu[k, ] <- .weighted_mean(params$data, w)
        self$sigma[[k]] <- .weighted_covariance(params$data, w, self$mu[k, ])
      }
      self$compute_density(params)
    }
  )
)

PoissonModel <- R6::R6Class(
  "PoissonModel",
  inherit = EmissionModel,
  public = list(
    lambda = NULL,

    initialize = function(params, Pi, post.pi) {
      self$lambda <- matrix(1, params$n_states, ncol(params$data))
      self$update(params, Pi, post.pi)
      .copy_init_fields(self, params$init$emission)
      self$compute_density(params)
    },

    compute_density = function(params) {
      K <- params$n_states
      P <- ncol(params$data)
      self$density <- matrix(0, params$n_obs, K)
      for (k in seq_len(K)) {
        log_d <- rep(0, params$n_obs)
        for (j in seq_len(P)) {
          log_d <- log_d + stats::dpois(
            params$data[, j], self$lambda[k, j], log = TRUE
          )
        }
        self$density[, k] <- exp(pmax(log_d, -700))
      }
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      for (k in seq_len(params$n_states)) {
        w <- post.pi[, k]
        sw <- max(sum(w), 1e-10)
        self$lambda[k, ] <- pmax(colSums(params$data * w) / sw, 1e-6)
      }
      self$compute_density(params)
    }
  )
)

ExponentialModel <- R6::R6Class(
  "ExponentialModel",
  inherit = EmissionModel,
  public = list(
    rate = NULL,

    initialize = function(params, Pi, post.pi) {
      self$rate <- matrix(1, params$n_states, ncol(params$data))
      self$update(params, Pi, post.pi)
      .copy_init_fields(self, params$init$emission)
      self$compute_density(params)
    },

    compute_density = function(params) {
      K <- params$n_states
      P <- ncol(params$data)
      self$density <- matrix(0, params$n_obs, K)
      for (k in seq_len(K)) {
        log_d <- rep(0, params$n_obs)
        for (j in seq_len(P)) {
          log_d <- log_d + stats::dexp(
            params$data[, j], self$rate[k, j], log = TRUE
          )
        }
        self$density[, k] <- exp(pmax(log_d, -700))
      }
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      for (k in seq_len(params$n_states)) {
        w <- post.pi[, k]
        mu <- .weighted_mean(params$data, w)
        self$rate[k, ] <- 1 / pmax(mu, 1e-6)
      }
      self$compute_density(params)
    }
  )
)

GammaModel <- R6::R6Class(
  "GammaModel",
  inherit = EmissionModel,
  public = list(
    shape = NULL,
    rate = NULL,

    initialize = function(params, Pi, post.pi) {
      K <- params$n_states
      P <- ncol(params$data)
      self$shape <- matrix(2, K, P)
      self$rate <- matrix(1, K, P)
      self$update(params, Pi, post.pi)
      .copy_init_fields(self, params$init$emission)
      self$compute_density(params)
    },

    compute_density = function(params) {
      K <- params$n_states
      P <- ncol(params$data)
      self$density <- matrix(0, params$n_obs, K)
      for (k in seq_len(K)) {
        log_d <- rep(0, params$n_obs)
        for (j in seq_len(P)) {
          log_d <- log_d + stats::dgamma(
            params$data[, j], self$shape[k, j], self$rate[k, j], log = TRUE
          )
        }
        self$density[, k] <- exp(pmax(log_d, -700))
      }
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      for (k in seq_len(params$n_states)) {
        w <- post.pi[, k]
        mu <- .weighted_mean(params$data, w)
        centered <- sweep(params$data, 2, mu, "-")
        va <- colSums(centered^2 * w) / max(sum(w), 1e-10)
        va <- pmax(va, 1e-6)
        self$shape[k, ] <- pmin(pmax(mu^2 / va, 0.05), 1000)
        self$rate[k, ] <- pmin(pmax(mu / va, 1e-6), 1000)
      }
      self$compute_density(params)
    }
  )
)

BetaModel <- R6::R6Class(
  "BetaModel",
  inherit = EmissionModel,
  public = list(
    shape1 = NULL,
    shape2 = NULL,

    initialize = function(params, Pi, post.pi) {
      K <- params$n_states
      P <- ncol(params$data)
      self$shape1 <- matrix(2, K, P)
      self$shape2 <- matrix(2, K, P)
      self$update(params, Pi, post.pi)
      .copy_init_fields(self, params$init$emission)
      self$compute_density(params)
    },

    compute_density = function(params) {
      K <- params$n_states
      P <- ncol(params$data)
      x <- pmin(pmax(params$data, 1e-8), 1 - 1e-8)
      self$density <- matrix(0, params$n_obs, K)
      for (k in seq_len(K)) {
        log_d <- rep(0, params$n_obs)
        for (j in seq_len(P)) {
          log_d <- log_d + stats::dbeta(
            x[, j], self$shape1[k, j], self$shape2[k, j], log = TRUE
          )
        }
        self$density[, k] <- exp(pmax(log_d, -700))
      }
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      x <- pmin(pmax(params$data, 1e-6), 1 - 1e-6)
      for (k in seq_len(params$n_states)) {
        w <- post.pi[, k]
        mu <- .weighted_mean(x, w)
        centered <- sweep(x, 2, mu, "-")
        va <- colSums(centered^2 * w) / max(sum(w), 1e-10)
        common <- mu * (1 - mu) / pmax(va, 1e-6) - 1
        common <- pmin(pmax(common, 0.1), 1000)
        self$shape1[k, ] <- pmax(mu * common, 0.05)
        self$shape2[k, ] <- pmax((1 - mu) * common, 0.05)
      }
      self$compute_density(params)
    }
  )
)

StudentTModel <- R6::R6Class(
  "StudentTModel",
  inherit = EmissionModel,
  public = list(
    mu = NULL,
    scale = NULL,
    df = NULL,

    initialize = function(params, Pi, post.pi) {
      K <- params$n_states
      P <- ncol(params$data)
      self$mu <- matrix(0, K, P)
      self$scale <- matrix(1, K, P)
      self$df <- rep(5, K)
      self$update(params, Pi, post.pi)
      .copy_init_fields(self, params$init$emission)
      self$compute_density(params)
    },

    compute_density = function(params) {
      K <- params$n_states
      P <- ncol(params$data)
      self$density <- matrix(0, params$n_obs, K)
      for (k in seq_len(K)) {
        log_d <- rep(0, params$n_obs)
        for (j in seq_len(P)) {
          z <- (params$data[, j] - self$mu[k, j]) / self$scale[k, j]
          log_d <- log_d + stats::dt(z, self$df[k], log = TRUE) -
            log(self$scale[k, j])
        }
        self$density[, k] <- exp(pmax(log_d, -700))
      }
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      for (k in seq_len(params$n_states)) {
        w <- post.pi[, k]
        self$mu[k, ] <- .weighted_mean(params$data, w)
        centered <- sweep(params$data, 2, self$mu[k, ], "-")
        va <- colSums(centered^2 * w) / max(sum(w), 1e-10)
        correction <- (self$df[k] - 2) / self$df[k]
        self$scale[k, ] <- sqrt(pmax(va * correction, 1e-6))
      }
      self$compute_density(params)
    }
  )
)
