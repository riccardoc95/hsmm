#' Fit a Hidden Semi-Markov Model (HSMM)
#'
#' Fits an HSMM (or HMM) from a matrix or data frame of observations while
#' optionally incorporating covariates for the transition and dwell-time
#' components and choosing among different emission families. The function
#' builds the duration, transition, and emission models via their respective
#' factories and runs the EM algorithm to estimate the parameters.
#'
#' @param data Matrix or data frame containing the observed sequence. Rows
#'   correspond to time points and columns to observed variables.
#' @param n_states Number of latent states of the model.
#' @param covariates_omega Optional covariates for the emission model
#'   \code{omega}. Must have the same number of rows as \code{data}.
#' @param covariates_q Optional covariates for the dwell-time model \code{q}.
#' @param semi Logical. If \code{TRUE} (default) the fitted model is an HSMM;
#'   if \code{FALSE} it reduces to a standard HMM.
#' @param max_dwell Maximum dwell time (in number of observations) allowed for
#'   each state in the semi-Markov formulation. If \code{NULL}, a default value
#'   of 1 is used.
#' @param model_type Name of the model specification to use for the emission
#'   component (for example \code{"torus"}, \code{"gaussian"}, etc.).
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance for the EM algorithm.
#' @param verbose Logical value indicating whether to print progress messages
#'   during training.
#' @param init Optional list with initial parameter values. If \code{NULL},
#'   starting values are determined automatically.
#' @param seed Numeric value used to set the random-number generator seed to
#'   ensure reproducibility.
#'
#' @return A list containing the estimated objects and the input parameters:
#' \itemize{
#'   \item \code{input_params}: parameters used during fitting.
#'   \item \code{emission.model}: fitted emission model object.
#'   \item \code{transition.model}: fitted transition model object.
#'   \item \code{duration.model}: fitted duration model object.
#' }
#'
#' @details The function initializes the maximum dwell times (if provided),
#' constructs the modular models through their factories, and launches the EM
#' algorithm to adapt the parameters to the supplied observations.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' obs <- data.frame(x = rnorm(100))
#' fit <- fit_hsmm(obs, n_states = 2, model_type = "torus",
#'                 max_iter = 10, verbose = FALSE)
#' }
#' @export
fit_hsmm <- function(
    # HMM parameters
  data,
  n_states,
  covariates_omega = NULL,
  covariates_q = NULL,

  # HSMM settings
  semi = TRUE,
  max_dwell = NULL,

  # Model Type
  model_type = "torus",

  # EM algorithm parameters
  max_iter = 100,
  tol = 1e-5,
  verbose = TRUE,
  init = NULL,
  seed = NULL) {

  # Collect function arguments into a parameter list
  params <- c(as.list(environment()))

  # Number of observations (rows in data)
  params$n_obs <- nrow(params$data)

  # ============================
  # Dwell-time and state indices
  # ============================
  # If user specifies a maximum dwell time, assign per-state dwell lengths.
  if(!is.null(params$max_dwell)){
    params$dwell_lengths <- rep(params$max_dwell, params$n_states)
    params$state_indices <- rep(1:params$n_states, params$dwell_lengths)
  } else {
    # Default: single-step dwell times (reduces to a standard HMM)
    params$max_dwell <- 1
    params$dwell_lengths <- rep(params$max_dwell, params$n_states)
    params$state_indices <- rep(1:params$n_states, params$dwell_lengths)
  }

  # ============================
  # Model construction (via factories)
  # ============================

  # Duration model
  duration.model.factory <- DurationModelFactory$new()
  duration.model <- duration.model.factory$create(params)

  # Transition model
  transition.model.factory <- TransitionModelFactory$new()
  transition.model <- transition.model.factory$create(params,
                                                      duration.model$p.array)

  # Emission model
  emission.model.factory <- EmissionModelFactory$new()
  emission.model <- emission.model.factory$create(params,
                                                  transition.model$Pi,
                                                  transition.model$post.pi)

  # ============================
  # EM Algorithm
  # ============================
  # Perform iterative parameter estimation
  llk <- em.algorithm(params,
                      emission.model, transition.model, duration.model)

  # ============================
  # Output
  # ============================
  # Return fitted model components and input configuration
  return(list(
    input_params = params,
    emission.model = emission.model,
    transition.model = transition.model,
    duration.model = duration.model,
    loglik = llk
  ))
}
