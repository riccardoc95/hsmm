#' Homogeneous Transition Model for HMM/HSMM
#'
#' Implements a *homogeneous* transition model, where transition probabilities
#' between hidden states remain constant over time (i.e., they do not depend
#' on external covariates or time-varying factors).
#'
#' @description
#' This model extends the base \code{TransitionModel} class. It builds a fixed
#' transition probability matrix \code{omega} and the corresponding augmented
#' transition tensor \code{gamma}, used for HSMM formulations that include
#' dwell-time modeling.
#'
#' @field omega  State-to-state transition probability matrix (K x K).
#' @field gamma  3D array of transition probabilities between augmented states (time × S × S).
#'
#' @details
#' - The model initializes \code{omega} as a valid stochastic matrix with no self-transitions.
#' - Transition probabilities within dwell times are set using the duration probability array (\code{p.array}).
#' - The final augmented transition tensor \code{gamma} is updated accordingly.
#'
#' @examples
#' params <- list(n_states = 3, max_dwell = 2, n_obs = 10)
#' p.array <- matrix(runif(3 * 2 * 10), nrow = 10, ncol = 6)
#' tm <- TransitionModel.Homogeneous$new(params, p.array)
#' tm$omega
#'
library(R6)

TransitionModel.Homogeneous <- R6Class(
  "TransitionModel.Homogeneous",
  inherit = TransitionModel,
  public = list(
    #' @description
    #' Initializes a homogeneous transition model with constant probabilities.
    #'
    #' @param params List of parameters containing:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states (K).}
    #'     \item{\code{n_obs}}{Number of time steps.}
    #'   }
    #' @param p.array Duration probability array (n × S) used to initialize \code{gamma}.
    initialize = function(params, p.array) {
      super$initialize(params)

      K <- params$n_states  # Number of hidden states

      # ---- Initialize transition matrix (omega) ----
      if (K == 1) {
        # Trivial case: single state (no transitions)
        self$omega <- matrix(1, 1, 1)

      } else if (K == 2) {
        # Simple deterministic transitions between 2 states
        self$omega <- matrix(c(0, 1,
                               1, 0), nrow = 2, byrow = TRUE)

      } else {
        # For K > 2, randomly initialize transition probabilities
        raw_omega <- matrix(runif(K * K), K, K)
        diag(raw_omega) <- 0  # No self-transitions

        # Normalize each row to sum to 1 (excluding diagonal)
        for (k in 1:K) {
          row_sum <- sum(raw_omega[k, -k])
          if (row_sum == 0) {
            # Handle potential numerical issue (avoid division by zero)
            raw_omega[k, -k] <- 1 / (K - 1)
          } else {
            raw_omega[k, -k] <- raw_omega[k, -k] / row_sum
          }
          raw_omega[k, k] <- 0
        }
        self$omega <- raw_omega
      }

      # Compute initial transition tensor gamma
      self$compute_gamma(params, p.array)
    },

    #' @description
    #' Computes the augmented transition probability tensor (\code{gamma})
    #' for all time steps and states.
    #'
    #' @param params List containing:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states (K).}
    #'     \item{\code{max_dwell}}{Maximum dwell time per state (M).}
    #'   }
    #' @param p.array Duration probability array, defining exit probabilities from dwell states.
    #' @return Invisibly returns \code{self}.
    compute_gamma = function(params, p.array) {
      K <- params$n_states
      M <- params$max_dwell

      self$gamma <- compute_gamma(self$omega,
                                  array(0, dim = c(1,1,1)),
                                  p.array, K, M, FALSE)


      # n     <- nrow(p.array)            # Number of observations
      # S_ext <- K * M                    # Total number of augmented states
      # dmap  <- rep(1:K, each = M)       # Mapping from augmented to original states
      #
      # # Initialize gamma as a 3D array: (time × from-state × to-state)
      # Gamma <- array(0, dim = c(n, S_ext, S_ext))
      #
      # # ---- Build gamma matrix for each time step ----
      # for (t in 1:n) {
      #   for (k in 1:K) {
      #     pos_k <- which(dmap == k)  # Augmented states belonging to main state k
      #
      #     # Within-state transitions (progressing through dwell time)
      #     if (length(pos_k) > 1) {
      #       for (m in 1:(length(pos_k) - 1)) {
      #         from_idx <- pos_k[m]
      #         to_idx   <- pos_k[m + 1]
      #         Gamma[t, from_idx, to_idx] <- 1 - p.array[t, from_idx]
      #       }
      #     }
      #
      #     # Transitions from last dwell state to other states
      #     last_sub <- tail(pos_k, 1)
      #     for (i in setdiff(1:K, k)) {
      #       pos_i <- which(dmap == i)
      #       first_i <- pos_i[1]
      #       Gamma[t, last_sub, first_i] <- p.array[t, last_sub] * self$omega[k, i]
      #     }
      #   }
      # }
      #
      # self$gamma <- Gamma
    },

    #' @description
    #' Updates the transition probability matrix (\code{omega})
    #' using posterior transition probabilities (\code{post.bi.pi.aug}).
    #'
    #' @param params List containing model configuration:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states.}
    #'     \item{\code{max_dwell}}{Maximum dwell time per state.}
    #'   }
    #' @param post.bi.pi.aug 3D array of posterior probabilities
    #' for transitions between augmented states.
    #' @return Invisibly returns \code{self}.
    compute_omega = function(params, post.bi.pi.aug) {
      K    <- params$n_states
      M    <- params$max_dwell
      dmap <- rep(1:K, each = M)  # Augmented-state mapping

      new_omega <- matrix(0, K, K)

      # ---- Aggregate posterior transitions ----
      for (k in 1:K) {
        from_states <- which(dmap == k)
        for (i in setdiff(1:K, k)) {
          to_states <- which(dmap == i)
          new_omega[k, i] <- sum(post.bi.pi.aug[, from_states, to_states, drop = FALSE])
        }

        # Normalize transition probabilities for each row
        row_sum <- sum(new_omega[k, -k])
        if (row_sum > 0) {
          new_omega[k, -k] <- new_omega[k, -k] / row_sum
        } else {
          # Default uniform transitions if no data support
          new_omega[k, -k] <- 1 / (K - 1)
        }
        new_omega[k, k] <- 0  # No self-transitions
      }

      self$omega <- new_omega
    }
  )
)
