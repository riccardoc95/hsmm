#' EM Step for Hidden Semi-Markov Model (HSMM)
#'
#' Performs one Expectation-Maximization (EM) step for a Hidden Semi-Markov Model (HSMM),
#' using emission, transition, and duration models provided by the user.
#'
#' @param params A list containing model parameters and settings, including:
#'   \itemize{
#'     \item{\code{n_states}}{Number of hidden states.}
#'     \item{\code{max_dwell}}{Maximum dwell time (duration) per state.}
#'     \item{\code{max_iter}}{Maximum number of EM iterations.}
#'     \item{\code{tol}}{Convergence tolerance for log-likelihood.}
#'     \item{\code{verbose}}{Logical, print progress if TRUE.}
#'   }
#' @param emission.model A list of functions defining the emission model,
#'   including:
#'   \itemize{
#'     \item{\code{density}}{Matrix of emission densities for each observation and state.}
#'     \item{\code{update}}{Function to update emission parameters.}
#'   }
#' @param transition.model A list defining the transition dynamics,
#'   with elements:
#'   \itemize{
#'     \item{\code{gamma}}{3D array of transition probabilities (time x from x to).}
#'     \item{\code{Pi}}{Initial state probability vector.}
#'     \item{\code{compute_omega, compute_gamma, update}}{Functions to update transition parameters.}
#'   }
#' @param duration.model A list of duration-related functions:
#'   \itemize{
#'     \item{\code{compute_beta, compute_p.array}}{Functions to update dwell-time distributions.}
#'     \item{\code{p.array}}{Array of duration probabilities.}
#'   }
#'
#' @return Numeric scalar — the log-likelihood (\code{l}) computed during this EM step.
#'
#' @details
#' This function implements the standard EM algorithm:
#' \itemize{
#'   \item \strong{E-step:} Compute expected sufficient statistics using the forward-backward algorithm.
#'   \item \strong{M-step:} Update model parameters based on posterior distributions.
#' }
#'
#' The helper function \code{backward.forward()} is used for the E-step.
#'
#' @seealso \code{\link{em.algorithm}}, \code{\link{backward.forward}}
em.step <- function(params, emission.model, transition.model, duration.model) {
  # =======================
  # E-STEP
  # =======================
  
  density <- emission.model$density               # Matrix of emission densities
  n_obs <- nrow(density)                          # Number of observations
  S <- params$n_states                            # Number of hidden states
  M <- params$max_dwell                           # Maximum dwell time per state
  SM <- S * M                                     # Total augmented state space size
  
  # Build emission probability matrix for augmented states
  fit <- matrix(0, nrow = n_obs, ncol = SM)
  for (k in 1:S) {
    pos <- ((k - 1) * M + 1):(k * M)              # Position indices for state k
    fit[, pos] <- matrix(rep(density[, k], M),
                         nrow = n_obs, ncol = M)  # Repeat densities for dwell times
  }
  
  # Replace invalid or zero values with machine epsilon to avoid numerical issues
  fit[!is.finite(fit) | fit <= 0] <- .Machine$double.eps
  
  # Run forward-backward algorithm
  bw <- backward.forward(
    fit                = fit,
    Gamma              = transition.model$gamma,
    Pi                 = transition.model$Pi,
    K                  = params$n_states,
    M                  = params$max_dwell
  )
  
  # =======================
  # M-STEP
  # =======================
  # Update duration and transition models using posterior probabilities
  
  duration.model$compute_beta(params, bw$post.bi.pi.aug)   # Update duration-related beta
  duration.model$compute_p.array(params)                   # Update duration probability array
  
  transition.model$compute_omega(params, bw$post.bi.pi.aug) # Update omega (transition intensity)
  transition.model$compute_gamma(params, duration.model$p.array) # Update transition probabilities
  transition.model$update(params, bw$post.pi.aug)          # Update transition model parameters
  
  # Collapse augmented posterior into standard state posterior
  n <- nrow(bw$post.pi.aug)
  post.pi.collapsed <- matrix(0, nrow = n, ncol = S)
  for (k in 1:S) {
    pos <- ((k - 1) * M + 1):(k * M)
    post.pi.collapsed[, k] <- rowSums(bw$post.pi.aug[, pos, drop = FALSE])
  }
  
  # Update emission model parameters using collapsed posterior
  emission.model$update(
    params,
    transition.model$Pi,          # Initial probabilities
    post.pi.collapsed             # Posterior probabilities per state
  )
  
  # Return log-likelihood from the forward-backward step
  return(bw$l)
}


#' EM Algorithm for Hidden Semi-Markov Models (HSMM)
#'
#' Performs the full Expectation-Maximization (EM) algorithm for an HSMM.
#' Iteratively calls \code{em.step()} until convergence or until the
#' maximum number of iterations is reached.
#'
#' @inheritParams em.step
#'
#' @return A numeric vector of log-likelihood values over all EM iterations.
#'
#' @details
#' The algorithm tracks the log-likelihood (\code{llk}) at each iteration and
#' terminates when the absolute difference between successive values is less than
#' \code{params$tol}. If \code{params$verbose = TRUE}, progress is printed at each iteration.
#'
#' @seealso \code{\link{em.step}}, \code{\link{backward.forward}}
em.algorithm <- function(params,
                         emission.model,
                         transition.model,
                         duration.model) {
  step <- 1
  llk <- numeric(params$max_iter)    # Vector to store log-likelihoods
  converged <- FALSE
  old.llk <- -Inf                   # Initialize previous log-likelihood
  
  # Main EM loop
  while (!converged && step <= params$max_iter) {
    
    # Perform one EM step
    llk[step] <- em.step(params,
                         emission.model,
                         transition.model,
                         duration.model)
    
    # Check for convergence (based on change in log-likelihood)
    if (step > 1) {
      dif <- abs(llk[step] - old.llk)
      if (!is.na(dif) && dif < params$tol) {
        converged <- TRUE
      }
    }
    
    # Optional progress output
    if (params$verbose) {
      cat("iteration", step,
          "; loglik =", llk[step],
          if (step > 1) paste0("; diff=", dif) else "",
          "\n")
    }
    
    # Update for next iteration
    old.llk <- llk[step]
    step <- step + 1
  }
  
  # Return only the log-likelihood values that were actually computed
  llk[1:(step - 1)]
}
