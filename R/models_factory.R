#' Duration, Transition, and Emission Model Factories
#'
#' These factory classes instantiate specific model subclasses used in
#' Hidden Semi-Markov Models (HSMMs) or Hidden Markov Models (HMMs)
#' depending on the user-supplied parameters.
#'
#' @details
#' The factory pattern allows the framework to flexibly construct different
#' duration, transition, and emission model types without changing core code.
#' Each factory inspects the contents of the \code{params} list to decide
#' which specific R6 model class to instantiate.
#'
#' @seealso
#' \itemize{
#'   \item \code{DurationModel.Geometric}, \code{DurationModel.Gompertz}
#'   \item \code{TransitionModel.Homogeneous}, \code{TransitionModel.TimeVarying}
#'   \item \code{TorusModel}, \code{GaussianModel}, \code{PoissonModel},
#'         \code{ExponentialModel}, \code{GammaModel}, \code{BetaModel},
#'         \code{StudentTModel}, \code{VarGaussianModel}
#' }
#'
#' @examples
#' \dontrun{
#' params <- list(n_states = 2, model_type = "gaussian")
#' em_factory <- EmissionModelFactory$new()
#' em_model <- em_factory$create(params, Pi = c(0.5, 0.5), post.pi = matrix(runif(20), ncol = 2))
#' }
#'
#' @name ModelFactories
NULL


#' Factory for Duration Models
#'
#' Creates and returns a duration model depending on whether covariates are provided.
#' If no covariates are present, a geometric duration model is used; otherwise, a
#' Gompertz-based duration model is chosen.
#'
DurationModelFactory <- R6Class(
  "DurationModelFactory",
  public = list(
    #' @description
    #' Constructs a duration model instance.
    #' @param params List of parameters including possible \code{covariates_q}.
    #' @return An R6 object inheriting from \code{DurationModel}.
    create = function(params) {
      # Determine model type based on the presence of covariates
      if (is.null(params$covariates_q) || ncol(params$covariates_q) == 0) {
        type <- "geometric"
      } else {
        type <- "gompertz"
      }
      
      # Instantiate appropriate duration model
      switch(
        type,
        "geometric" = DurationModel.Geometric$new(params),
        "gompertz"  = DurationModel.Gompertz$new(params),
        stop("Unknown duration model type: ", type)
      )
    }
  )
)


#' Factory for Transition Models
#'
#' Creates and returns a transition model depending on the presence of covariates.
#' Without covariates, a homogeneous transition model is used; with covariates,
#' a time-varying transition model is constructed.
#'
TransitionModelFactory <- R6Class(
  "TransitionModelFactory",
  public = list(
    #' @description
    #' Constructs a transition model instance.
    #' @param params List of parameters, including \code{covariates_omega}.
    #' @param p.array Duration probability array (used by the transition model).
    #' @return An R6 object inheriting from \code{TransitionModel}.
    create = function(params, p.array) {
      # Decide model type based on covariates
      if (is.null(params$covariates_omega) || ncol(params$covariates_omega) == 0) {
        type <- "homogeneous"
      } else {
        type <- "timevarying"
      }
      
      # Instantiate appropriate transition model
      switch(
        type,
        "homogeneous" = TransitionModel.Homogeneous$new(params, p.array),
        "timevarying" = TransitionModel.TimeVarying$new(params, p.array),
        stop("Unknown transition model type: ", type)
      )
    }
  )
)


#' Factory for Emission Models
#'
#' Creates and returns an emission model instance based on the specified
#' \code{model_type} in the \code{params} list. Supports a variety of
#' emission families such as Gaussian, Poisson, Gamma, Beta, Torus, etc.
#'
EmissionModelFactory <- R6Class(
  "EmissionModelFactory",
  public = list(
    #' @description
    #' Constructs an emission model instance.
    #' @param params List of model parameters including \code{model_type}.
    #' @param Pi Numeric vector of initial state probabilities.
    #' @param post.pi Matrix of posterior probabilities per state.
    #' @return An R6 object inheriting from \code{EmissionModel}.
    create = function(params, Pi, post.pi) {
      # Convert model type to lowercase for safety
      type <- tolower(params$model_type)
      
      # Instantiate correct emission model based on model_type
      switch(
        type,
        "torus"      = TorusModel$new(params, Pi, post.pi),
        "gaussian"   = GaussianModel$new(params, Pi, post.pi),
        "poisson"    = PoissonModel$new(params, Pi, post.pi),
        "exponential"= ExponentialModel$new(params, Pi, post.pi),
        "gamma"      = GammaModel$new(params, Pi, post.pi),
        "beta"       = BetaModel$new(params, Pi, post.pi),
        "tstudent"   = StudentTModel$new(params, Pi, post.pi),
        "vargaussian"= VarGaussianModel$new(params, Pi, post.pi),
        stop("Unknown emission model type: ", type)
      )
    }
  )
)
