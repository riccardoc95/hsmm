# Simple state-specific Gaussian vector autoregression.
VarGaussianModel <- R6::R6Class(
  "VarGaussianModel",
  inherit = EmissionModel,
  public = list(
    ar = NULL,
    lambda = NULL,
    coef = NULL,
    sigma = NULL,
    rows = NULL,

    initialize = function(params, Pi, post.pi) {
      self$ar <- params$ar
      self$lambda <- params$lambda
      self$coef <- vector("list", params$n_states)
      self$sigma <- vector("list", params$n_states)
      self$update(params, Pi, post.pi)
      .copy_init_fields(self, params$init$emission)
      self$compute_density(params)
    },

    compute_density = function(params) {
      lagged <- .make_lagged_data(params$data, self$ar)
      self$rows <- lagged$rows
      self$density <- matrix(1, params$n_obs, params$n_states)

      for (k in seq_len(params$n_states)) {
        mu <- lagged$X %*% self$coef[[k]]
        residuals <- lagged$Y - mu
        zero <- rep(0, ncol(params$data))
        self$density[self$rows, k] <- .dmvnorm_basic(
          residuals, zero, self$sigma[[k]]
        )
      }
      self$density <- pmax(self$density, .Machine$double.xmin)
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      lagged <- .make_lagged_data(params$data, self$ar)
      for (k in seq_len(params$n_states)) {
        ans <- .weighted_var_fit(
          lagged$X, lagged$Y, post.pi[lagged$rows, k], self$lambda
        )
        self$coef[[k]] <- ans$coef
        self$sigma[[k]] <- ans$sigma
      }
      self$compute_density(params)
    }
  )
)
