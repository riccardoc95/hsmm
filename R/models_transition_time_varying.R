# Conditional transition probabilities described by multinomial logits.
TransitionModel.TimeVarying <- R6::R6Class(
  "TransitionModel.TimeVarying",
  inherit = TransitionModel,
  public = list(
    initialize = function(params, p.array) {
      super$initialize(params)
      K <- params$n_states
      X <- cbind(1, params$covariates_omega[seq_len(params$n_obs - 1), , drop = FALSE])

      self$coefficients <- vector("list", K)
      self$destinations <- vector("list", K)
      for (k in seq_len(K)) {
        dest <- setdiff(seq_len(K), k)
        self$destinations[[k]] <- dest
        self$coefficients[[k]] <- matrix(0, ncol(X), max(length(dest) - 1, 0))
      }

      initial <- params$init$transition_coefficients
      if (is.null(initial) && !is.null(params$init$transition)) {
        initial <- params$init$transition$coefficients
      }
      if (!is.null(initial)) self$coefficients <- initial

      self$update_omega(params)
      self$compute_gamma(params, p.array)
    },

    update_omega = function(params) {
      K <- params$n_states
      Tm1 <- params$n_obs - 1
      X <- cbind(1, params$covariates_omega[seq_len(Tm1), , drop = FALSE])
      self$omega_time <- array(0, dim = c(Tm1, K, K))

      if (K == 1) {
        self$omega_time[, 1, 1] <- 1
        self$omega <- matrix(1, 1, 1)
        return(invisible(self))
      }

      for (k in seq_len(K)) {
        dest <- self$destinations[[k]]
        pr <- .predict_softmax_reference(
          X, self$coefficients[[k]], length(dest)
        )
        for (j in seq_along(dest)) {
          self$omega_time[, k, dest[j]] <- pr[, j]
        }
      }

      self$omega <- apply(self$omega_time, c(2, 3), mean)
      diag(self$omega) <- 0
      invisible(self)
    },

    compute_gamma = function(params, p.array) {
      self$gamma <- compute_gamma(
        matrix(0, params$n_states, params$n_states),
        self$omega_time,
        p.array,
        params$n_states,
        params$max_dwell,
        TRUE
      )
      invisible(self)
    },

    compute_omega = function(params, post.bi.pi.aug) {
      K <- params$n_states
      M <- params$max_dwell
      Tm1 <- params$n_obs - 1
      if (K <= 2) return(invisible(self))

      X <- cbind(1, params$covariates_omega[seq_len(Tm1), , drop = FALSE])
      for (k in seq_len(K)) {
        Y <- .transition_counts(post.bi.pi.aug, K, M, k)
        self$coefficients[[k]] <- .fit_transition_softmax(
          X, Y, self$coefficients[[k]]
        )
      }
      self$update_omega(params)
      invisible(self)
    }
  )
)
