# Base class for transitions between the latent states.
TransitionModel <- R6::R6Class(
  "TransitionModel",
  public = list(
    gamma = NULL,
    omega = NULL,
    omega_time = NULL,
    coefficients = NULL,
    destinations = NULL,
    Pi = NULL,
    post.pi = NULL,

    initialize = function(params) {
      K <- params$n_states
      self$post.pi <- params$initial_post
      self$Pi <- .normalize_probability(self$post.pi[1, ])

      if (!is.null(params$init$Pi)) {
        if (length(params$init$Pi) != K) stop("init$Pi must have length n_states")
        self$Pi <- .normalize_probability(params$init$Pi)
      }
    },

    compute_gamma = function(params, p.array) {
      invisible(self)
    },

    compute_omega = function(params, post.bi.pi.aug) {
      invisible(self)
    },

    update = function(params, post.pi.aug) {
      self$post.pi <- .collapse_posterior(
        post.pi.aug, params$n_states, params$max_dwell
      )
      self$Pi <- .normalize_probability(self$post.pi[1, ])
      invisible(self)
    }
  )
)
