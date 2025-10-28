#' Base Transition Model for HMM/HSMM
#'
#' Defines the base structure for the transition component in Hidden Markov
#' Models (HMMs) and Hidden Semi-Markov Models (HSMMs).  
#' This class manages the transition probability matrices (`gamma`),
#' transition intensities (`omega`), and initial state probabilities (`Pi`).
#'
#' @field gamma  3D array of transition probabilities between augmented states.
#' @field omega  Transition intensity parameters (may vary by model subclass).
#' @field Pi     Numeric vector of initial state probabilities.
#' @field post.pi Matrix of posterior probabilities per state (updated during EM).
#'
#' @details
#' This base class provides general initialization and update logic for state
#' probabilities. Specific model subclasses (e.g., homogeneous or time-varying)
#' should override \code{compute_gamma()}, \code{compute_omega()}, and possibly
#' the \code{update()} method to implement their custom transition structures.
#'
#' @examples
#' params <- list(n_states = 3, n_obs = 10, state_indices = rep(1:3, each = 1))
#' tm <- TransitionModel$new(params)
#' tm$Pi  # Initial state probabilities
#'
library(R6)

TransitionModel <- R6Class(
  "TransitionModel",
  public = list(
    gamma   = NULL,  # Transition probability tensor (n-1 x S x S)
    omega   = NULL,  # Transition intensity or auxiliary parameters
    Pi      = NULL,  # Initial state distribution
    post.pi = NULL,  # Posterior state probabilities
    
    #' @description
    #' Initializes the transition model with uniform probabilities.
    #' 
    #' @param params A list containing model configuration, including:
    #' \itemize{
    #'   \item{\code{n_states}}{Number of hidden states.}
    #'   \item{\code{n_obs}}{Number of time steps (observations).}
    #' }
    initialize = function(params) {
      K <- params$n_states   # Number of states
      n <- params$n_obs      # Number of observations
      
      # Initialize uniform state probabilities
      self$Pi <- rep(1 / K, K)
      self$post.pi <- matrix(1 / K, nrow = n, ncol = K)
    },
    
    #' @description
    #' Computes the transition probability array (\code{gamma}).
    #' Placeholder method — intended to be overridden in subclasses.
    #' 
    #' @param params Model parameters.
    #' @param p.array Duration probability array (used in HSMMs).
    #' @return Invisibly returns \code{self}.
    compute_gamma = function(params, p.array) {
      invisible(self)
    },
    
    #' @description
    #' Computes the transition intensity or auxiliary parameter matrix (\code{omega}).
    #' Placeholder method — intended for subclasses.
    #' 
    #' @param params Model parameters.
    #' @param post.bi.pi.aug Posterior probabilities of state transitions.
    #' @return Invisibly returns \code{self}.
    compute_omega = function(params, post.bi.pi.aug) {
      invisible(self)
    },
    
    #' @description
    #' Updates the posterior state probabilities and initial state distribution
    #' based on augmented posterior probabilities.
    #' 
    #' @param params List containing model configuration (\code{n_states}, \code{n_obs}, \code{state_indices}).
    #' @param post.pi.aug Matrix of posterior probabilities for augmented states.
    update = function(params, post.pi.aug) {
      K  <- params$n_states         # Number of latent states
      n  <- params$n_obs            # Number of time steps
      si <- params$state_indices    # Mapping from augmented states to main states
      
      # Aggregate augmented posterior probabilities into per-state posterior
      self$post.pi <- vapply(
        1:n,
        function(t) {
          tapply(post.pi.aug[t, ], si, sum)
        },
        numeric(K)
      )
      self$post.pi <- t(self$post.pi)
      
      # Update initial state probabilities
      self$Pi <- self$post.pi[1, ]
    }
  )
)
