# Independent von Mises margins for angular observations.
TorusModel <- R6::R6Class(
  "TorusModel",
  inherit = EmissionModel,
  public = list(
    mu = NULL,
    kappa = NULL,

    initialize = function(params, Pi, post.pi) {
      K <- params$n_states
      P <- ncol(params$data)
      self$mu <- matrix(0, K, P)
      self$kappa <- matrix(1, K, P)
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
          kap <- self$kappa[k, j]
          scaled_bessel <- besselI(kap, 0, expon.scaled = TRUE)
          log_norm <- log(2 * pi) + log(scaled_bessel) + kap
          log_d <- log_d + kap * cos(params$data[, j] - self$mu[k, j]) -
            log_norm
        }
        self$density[, k] <- exp(pmax(log_d, -700))
      }
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      for (k in seq_len(params$n_states)) {
        w <- post.pi[, k]
        sw <- max(sum(w), 1e-10)
        for (j in seq_len(ncol(params$data))) {
          cbar <- sum(w * cos(params$data[, j])) / sw
          sbar <- sum(w * sin(params$data[, j])) / sw
          self$mu[k, j] <- atan2(sbar, cbar) %% (2 * pi)
          self$kappa[k, j] <- .kappa_from_resultant(sqrt(cbar^2 + sbar^2))
        }
      }
      self$compute_density(params)
    }
  )
)
