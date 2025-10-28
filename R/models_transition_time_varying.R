#' Time-Varying Transition Model for HMM/HSMM
#'
#' Implements a *time-varying* transition model, where transition probabilities
#' between hidden states are allowed to change across time steps.
#'
#' @description
#' This class extends \code{TransitionModel} and generalizes the homogeneous case by
#' allowing a separate transition probability matrix (\code{omega[t,,]}) for each time step.
#' This makes it possible to model systems where state transitions evolve dynamically
#' over time or depend implicitly on external covariates.
#'
#' @field omega  3D array of transition probabilities (\code{[time × from × to]}).
#' @field gamma  3D array of transition probabilities between augmented states.
#'
#' @details
#' - When \code{n_states > 2}, each row of \code{omega[t,,]} is normalized to form
#'   a valid probability distribution without self-transitions.
#' - Within-state transitions (due to dwell time) are computed from the
#'   duration probability array (\code{p.array}).
#' - Cross-state transitions depend on the corresponding time-specific
#'   matrix \code{omega[t,,]}.
#'
#' @examples
#' params <- list(n_states = 3, max_dwell = 2, n_obs = 10)
#' p.array <- matrix(runif(10 * 6), nrow = 10)
#' tm <- TransitionModel.TimeVarying$new(params, p.array)
#' dim(tm$omega)
#'
library(R6)

TransitionModel.TimeVarying <- R6Class(
  "TransitionModel.TimeVarying",
  inherit = TransitionModel,
  public = list(
    #' @description
    #' Initializes a time-varying transition model.
    #' Creates a separate transition matrix for each time step (\code{omega[t,,]}).
    #' 
    #' @param params List containing:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states.}
    #'     \item{\code{n_obs}}{Number of time steps.}
    #'   }
    #' @param p.array Duration probability array used to initialize \code{gamma}.
    initialize = function(params, p.array) {
      super$initialize(params)
      
      K <- params$n_states     # Number of hidden states
      Tm1 <- params$n_obs - 1  # Number of transitions (time steps minus one)
      
      # Initialize 3D transition array (time × from × to)
      self$omega <- array(0, dim = c(Tm1, K, K))
      
      # ---- Initialize time-specific omega matrices ----
      if (K == 1) {
        # Trivial case: single-state system
        for (t in 1:Tm1) {
          self$omega[t, , ] <- matrix(1, 1, 1)
        }
        
      } else if (K == 2) {
        # Deterministic flipping between 2 states at all time steps
        for (t in 1:Tm1) {
          self$omega[t, , ] <- matrix(c(0, 1,
                                        1, 0), nrow = 2, byrow = TRUE)
        }
        
      } else {
        # Random initialization for each time step (K > 2)
        for (t in 1:Tm1) {
          raw_omega_t <- matrix(runif(K * K), K, K)
          diag(raw_omega_t) <- 0  # No self-transitions
          
          # Normalize each row
          for (k in 1:K) {
            row_sum <- sum(raw_omega_t[k, -k])
            if (row_sum == 0) {
              raw_omega_t[k, -k] <- 1 / (K - 1)
            } else {
              raw_omega_t[k, -k] <- raw_omega_t[k, -k] / row_sum
            }
            raw_omega_t[k, k] <- 0
          }
          self$omega[t, , ] <- raw_omega_t
        }
      }
      
      # Compute initial transition tensor gamma
      self$compute_gamma(params, p.array)
    },
    
    #' @description
    #' Computes the augmented transition tensor (\code{gamma}) for all time steps,
    #' taking into account time-varying transition probabilities.
    #' 
    #' @param params List of parameters containing:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states.}
    #'     \item{\code{max_dwell}}{Maximum dwell time per state.}
    #'   }
    #' @param p.array Duration probability array.
    #' @return Invisibly returns \code{self}.
    compute_gamma = function(params, p.array) {
      K <- params$n_states
      M <- params$max_dwell
      
      n     <- nrow(p.array)        # Number of time steps
      S_ext <- K * M                # Number of augmented states
      dmap  <- rep(1:K, each = M)   # Map from augmented to main states
      
      # Initialize gamma tensor (time × from × to)
      Gamma <- array(0, dim = c(n, S_ext, S_ext))
      
      # ---- Build augmented transition matrix for each time step ----
      for (t in 1:n) {
        for (k in 1:K) {
          pos_k <- which(dmap == k)  # Augmented states for state k
          
          # Within-state dwell transitions
          if (length(pos_k) > 1) {
            for (m in 1:(length(pos_k) - 1)) {
              from_idx <- pos_k[m]
              to_idx   <- pos_k[m + 1]
              Gamma[t, from_idx, to_idx] <- 1 - p.array[t, from_idx]
            }
          }
          
          # Between-state transitions (time-varying)
          last_sub <- tail(pos_k, 1)
          for (i in setdiff(1:K, k)) {
            pos_i   <- which(dmap == i)
            first_i <- pos_i[1]
            
            # Use omega for time t (or last available if t > Tm1)
            Gamma[t, last_sub, first_i] <- p.array[t, last_sub] *
              self$omega[min(t, dim(self$omega)[1]), k, i]
          }
        }
      }
      
      self$gamma <- Gamma
    },
    
    #' @description
    #' Updates the time-varying transition matrices (\code{omega[t,,]})
    #' using posterior transition probabilities.
    #' 
    #' @param params List of model parameters:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states.}
    #'     \item{\code{max_dwell}}{Maximum dwell time per state.}
    #'     \item{\code{n_obs}}{Number of time steps.}
    #'   }
    #' @param post.bi.pi.aug 3D array of posterior transition probabilities
    #' between augmented states.
    #' @return Invisibly returns \code{self}.
    compute_omega = function(params, post.bi.pi.aug) {
      K    <- params$n_states
      M    <- params$max_dwell
      Tm1  <- params$n_obs - 1
      dmap <- rep(1:K, each = M)
      
      # Initialize updated omega tensor
      new_omega <- array(0, dim = c(Tm1, K, K))
      
      # ---- Update omega for each time step ----
      for (t in 1:Tm1) {
        for (k in 1:K) {
          from_states <- which(dmap == k)
          for (i in setdiff(1:K, k)) {
            to_states <- which(dmap == i)
            new_omega[t, k, i] <- sum(post.bi.pi.aug[t, from_states, to_states, drop = FALSE])
          }
          
          # Normalize each row to form a valid probability distribution
          row_sum <- sum(new_omega[t, k, -k])
          if (row_sum > 0) {
            new_omega[t, k, -k] <- new_omega[t, k, -k] / row_sum
          } else {
            new_omega[t, k, -k] <- 1 / (K - 1)
          }
          new_omega[t, k, k] <- 0
        }
      }
      
      self$omega <- new_omega
    }
  )
)
