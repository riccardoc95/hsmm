#' Base Emission Model for HMM/HSMM
#'
#' Defines the abstract interface for emission (observation) models used
#' in Hidden Markov Models (HMMs) and Hidden Semi-Markov Models (HSMMs).
#'
#' @description
#' The \code{EmissionModel} class provides the base structure and generic methods
#' for modeling the conditional distribution of observations given hidden states.
#' 
#' Specific subclasses (e.g., \code{GaussianModel}, \code{PoissonModel},
#' \code{GammaModel}, \code{TorusModel}) must override:
#' \itemize{
#'   \item \code{compute_density()} — to compute the probability density
#'         or mass function for each observation and state.
#'   \item \code{update()} — to re-estimate emission parameters during the M-step
#'         of the EM algorithm based on posterior state probabilities.
#' }
#'
#' @field density A numeric matrix or array storing the emission densities
#'   for each observation (rows) and each hidden state (columns).
#'
#' @examples
#' # Example of subclassing the EmissionModel
#' SimpleGaussianModel <- R6::R6Class(
#'   "SimpleGaussianModel",
#'   inherit = EmissionModel,
#'   public = list(
#'     mu = NULL,
#'     sigma = NULL,
#'     compute_density = function(params, p.array) {
#'       self$density <- dnorm(params$data, mean = self$mu, sd = self$sigma)
#'     }
#'   )
#' )
#'
library(R6)

EmissionModel <- R6Class(
  "EmissionModel",
  public = list(
    #' @field density Numeric matrix or array of emission densities
    #'   (observations × states).
    density = NULL,
    
    #' @description
    #' Base initialization method for the emission model.
    #' 
    #' @param params List of model parameters (data, number of states, etc.).
    #' @param Pi Initial state probability vector.
    #' @param post.pi Posterior probabilities from the previous EM iteration.
    #' @return Invisibly returns \code{self}.
    initialize = function(params, Pi, post.pi) {
      # This base class does not initialize any parameters by default.
      # Subclasses should set up model-specific parameters (e.g., mean, variance).
      invisible(self)
    },
    
    #' @description
    #' Placeholder for computing emission densities.
    #' Subclasses must override this method to implement
    #' the actual density computation logic.
    #' 
    #' @param params Model parameter list.
    #' @param p.array Optional additional parameter (for compatibility).
    #' @return Invisibly returns \code{self}.
    compute_density = function(params, p.array) {
      # Should compute and store emission densities (per observation × state)
      invisible(self)
    },
    
    #' @description
    #' Placeholder for updating emission model parameters during the M-step.
    #' Subclasses must override this method to implement
    #' maximum-likelihood parameter updates.
    #' 
    #' @param params Model parameter list.
    #' @param Pi Current initial state probability vector.
    #' @param post.pi Posterior state probabilities.
    #' @return Invisibly returns \code{self}.
    update = function(params, Pi, post.pi) {
      # Should re-estimate emission parameters based on posterior weights.
      invisible(self)
    }
  )
)
