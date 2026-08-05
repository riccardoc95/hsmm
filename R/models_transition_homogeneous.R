# Transition probabilities that are constant through time.
TransitionModel.Homogeneous <- R6::R6Class(
  "TransitionModel.Homogeneous",
  inherit = TransitionModel,
  public = list(
    initialize = function(params, p.array) {
      super$initialize(params)
      K <- params$n_states

      if (K == 1) {
        self$omega <- matrix(1, 1, 1)
      } else if (K == 2) {
        self$omega <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
      } else {
        self$omega <- matrix(1 / (K - 1), K, K)
        diag(self$omega) <- 0
      }

      if (!is.null(params$init$omega)) self$omega <- params$init$omega
      if (!is.null(params$init$transition$omega)) {
        self$omega <- params$init$transition$omega
      }
      self$omega <- .normalize_omega(self$omega)
      self$compute_gamma(params, p.array)
    },

    compute_gamma = function(params, p.array) {
      self$gamma <- compute_gamma(
        self$omega,
        array(0, dim = c(1, 1, 1)),
        p.array,
        params$n_states,
        params$max_dwell,
        FALSE
      )
      invisible(self)
    },

    compute_omega = function(params, post.bi.pi.aug) {
      K <- params$n_states
      M <- params$max_dwell
      if (K <= 2) return(invisible(self))

      new_omega <- matrix(0, K, K)
      for (k in seq_len(K)) {
        dest <- setdiff(seq_len(K), k)
        counts <- .transition_counts(post.bi.pi.aug, K, M, k)
        totals <- colSums(counts)
        if (sum(totals) > 1e-12) {
          new_omega[k, dest] <- totals / sum(totals)
        } else {
          new_omega[k, dest] <- self$omega[k, dest]
        }
      }
      diag(new_omega) <- 0
      self$omega <- new_omega
      invisible(self)
    }
  )
)
