#' DurationModel R6 Class
#'
#' Represents the duration (or dwell-time) component of a Hidden Semi-Markov Model (HSMM).
#' This class defines the duration distribution parameters, including intercepts,
#' time effects, and optional covariate effects. It also provides methods for computing
#' the duration probability array and updating duration parameters during the EM algorithm.
#'
#' @details
#' The duration model defines the probability of remaining in a given state for
#' a certain number of time steps (dwell time). The underlying hazard function is
#' based on a **Gompertz distribution** using a complementary log-log link (\code{cloglog}).
#'
#' @field p.array Array storing computed duration probabilities for each state and dwell time.
#' @field beta.intercepts Numeric vector of intercept parameters for each state.
#' @field beta.time_effects Numeric vector of duration-dependent slope parameters.
#' @field beta.covariate Optional coefficients for duration covariates.
#' @field d.state Numeric vector indexing which state each duration corresponds to.
#' @field d.time Numeric vector of dwell-time values (1, 2, ..., M).
#' @field d.covariates Optional matrix of covariate values for dwell-time modeling.
#' @field link Character string indicating the link function used (default: "cloglog").
#'
#' @examples
#' params <- list(n_states = 2, max_dwell = 5, semi = TRUE)
#' dm <- DurationModel$new(params)
#' dm$beta.intercepts
#'
library(R6)

DurationModel <- R6Class(
  "DurationModel",
  public = list(
    # Public fields (model parameters)
    p.array = NULL,           # Duration probability array
    beta.intercepts = NULL,   # Intercept coefficients per state
    beta.time_effects = NULL, # Time-dependent coefficients per state
    beta.covariate = NULL,    # Optional covariate coefficients (unused here)
    
    d.state = NULL,           # State index associated with each dwell-time position
    d.time = NULL,            # Time indices (1:M) for each state
    d.covariates = NULL,      # Optional duration covariates
    
    link = NULL,              # Link function name (default "cloglog")
    
    #' @description
    #' Initializes a new DurationModel object.
    #' 
    #' @param params A list containing:
    #' \itemize{
    #'   \item{\code{n_states}}{Number of hidden states.}
    #'   \item{\code{max_dwell}}{Maximum dwell time per state.}
    #'   \item{\code{semi}}{Logical, if TRUE (default) use semi-Markov structure.}
    #' }
    initialize = function(params) {
      K <- params$n_states     # Number of hidden states
      M <- params$max_dwell    # Maximum dwell time per state
      
      # Random initialization of intercepts for each state
      self$beta.intercepts <- runif(K, -3, -0.5)
      
      # State index vector (repeated for each dwell-time level)
      self$d.state <- rep(1:K, each = M)
      
      if (params$semi) {
        # Semi-Markov case: include time effect on hazard
        self$beta.time_effects <- runif(K, 0.01, 0.3)
        # Time index vector (1:M) repeated for each state
        self$d.time <- rep(1:M, times = K)
      } else {
        # Standard HMM case: no time-dependent effect
        self$beta.time_effects <- rep(0, K)
        self$d.time <- rep(0, K * M)
      }
      
      # Use complementary log-log link for hazard modeling
      self$link <- "cloglog"
    },
    
    #' @description
    #' Gompertz hazard base function.
    #' Computes the hazard value for a given time index using a Gompertz-type
    #' parameterization: \eqn{h(i) = exp(alpha + beta * i)}.
    #'
    #' @param i Integer or numeric vector representing time indices (1:M).
    #' @param alpha Numeric intercept parameter.
    #' @param betai Numeric slope parameter for time effect.
    #' @return Numeric vector of hazard values.
    gomp = function(i, alpha, betai) {
      exp(pmin(alpha + betai * i, 700))  # Clamp exponent to 700 for numerical stability
    },
    
    #' @description
    #' Computes the duration probability array (\code{p.array}).
    #' Placeholder method — should be implemented in subclass models.
    #' 
    #' @param params Model parameter list.
    #' @return Invisibly returns \code{self} for method chaining.
    compute_p.array = function(params) {
      invisible(self)
    },
    
    #' @description
    #' Updates the duration model coefficients (\eqn{\beta}) based on
    #' posterior transition probabilities.
    #' Placeholder method — to be customized by specific duration models.
    #'
    #' @param params Model parameter list.
    #' @param post.bi.pi.aug 3D array of pairwise posterior transition probabilities.
    #' @return Invisibly returns \code{self} for method chaining.
    compute_beta = function(params, post.bi.pi.aug) {
      invisible(self)
    }
  )
)
