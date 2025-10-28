#' Torus Emission Model for HMM/HSMM
#'
#' Represents a bivariate angular (torus) emission model used for circular data,
#' where each observation lies on the product of two unit circles (e.g. angles).
#' The density is computed using a parametric torus distribution controlled by
#' parameters \eqn{\mu_1, \mu_2, \kappa_1, \kappa_2, \rho}.
#'
#' @field tau.mu1 Numeric vector of mean directions for the first angle (\eqn{\mu_1}).
#' @field tau.mu2 Numeric vector of mean directions for the second angle (\eqn{\mu_2}).
#' @field tau.k1 Numeric vector of concentration parameters for the first angle (\eqn{\kappa_1}).
#' @field tau.k2 Numeric vector of concentration parameters for the second angle (\eqn{\kappa_2}).
#' @field tau.rho Numeric vector of correlation parameters (\eqn{\rho}) between the two angles.
#' @field density Matrix of computed emission densities for each observation and state.
#'
#' @details
#' This model is suited for **bivariate circular data**, such as directional
#' statistics (e.g., wind direction pairs or phase angles).  
#' Each state's emission is parameterized using the torus distribution, where
#' the joint density is calculated from trigonometric components.
#'
#' @examples
#' \dontrun{
#' # Example with dummy circular data (angles in radians)
#' params <- list(
#'   n_states = 2,
#'   data = cbind(runif(100, 0, 2*pi), runif(100, 0, 2*pi))
#' )
#' Pi <- rep(1/2, 2)
#' post.pi <- matrix(runif(200), ncol = 2)
#' model <- TorusModel$new(params, Pi, post.pi)
#' }
#'
library(R6)

TorusModel <- R6Class(
  "TorusModel",
  inherit = EmissionModel,
  
  public = list(
    tau.mu1 = NULL,   # Mean direction parameter for first angle
    tau.mu2 = NULL,   # Mean direction parameter for second angle
    tau.k1 = NULL,    # Concentration parameter for first angle
    tau.k2 = NULL,    # Concentration parameter for second angle
    tau.rho = NULL,   # Correlation parameter between the two angles
    density = NULL,   # Matrix of emission densities
    
    #' @description
    #' Initializes the torus emission model parameters.
    #' Parameters are estimated from weighted trigonometric moments of the data.
    #' 
    #' @param params List containing model parameters, including \code{data} and \code{n_states}.
    #' @param Pi Initial state probabilities (optional).
    #' @param post.pi Posterior state probabilities (optional).
    initialize = function(params, Pi = NULL, post.pi = NULL) {
      super$initialize(params, Pi, post.pi)
      
      # Initialize parameters (tau.mu*, tau.k*, tau.rho) via circular moments
      theta <- matrix(0, nrow = params$n_states, ncol = 5)
      for (k in 1:params$n_states) {
        trig.mom <- colSums(post.pi[, k] * (cos(params$data) + 1i * sin(params$data))) / sum(post.pi[, k])
        theta[k, ] <- c(Arg(trig.mom[1]), Arg(trig.mom[2]), Mod(trig.mom[1]), Mod(trig.mom[2]), 0)
      }
      self$tau.mu1 <- theta[, 1]
      self$tau.mu2 <- theta[, 2]
      self$tau.k1 <- theta[, 3]
      self$tau.k2 <- theta[, 4]
      self$tau.rho <- theta[, 5]
      
      # Compute initial emission densities
      self$compute_density(params, Pi)
    },
    
    #' @description
    #' Computes the emission density for all observations and states.
    #' @param params List of model parameters (including data and number of states).
    #' @param Pi Optional vector of state probabilities.
    compute_density = function(params, Pi) {
      data <- params$data
      n_obs <- nrow(data)
      n_states <- self$n_states
      
      dens <- array(dim = c(n_obs, params$n_states))
      for (k in 1:params$n_states) {
        tau.par <- c(self$tau.mu1[k], self$tau.mu2[k],
                     self$tau.k1[k], self$tau.k2[k], self$tau.rho[k])
        dens[, k] <- self$density_function(data[, 1], data[, 2], Pi, tau.par)
      }
      
      self$density <- dens
      invisible(self)
    },
    
    #' @description
    #' Torus density function.
    #' Computes the joint circular density for given angles (\eqn{\theta_1}, \eqn{\theta_2})
    #' and parameters \eqn{(\mu_1, \mu_2, \kappa_1, \kappa_2, \rho)}.
    #' 
    #' @param theta1 Vector of first angle values.
    #' @param theta2 Vector of second angle values.
    #' @param Pi Unused here (for compatibility).
    #' @param par Numeric vector of length 5 containing the torus parameters.
    #' @return Numeric vector of density values.
    density_function = function(theta1, theta2, Pi, par = NULL) {
      tau.mu1 <- par[1]; tau.mu2 <- par[2]
      tau.k1 <- par[3]; tau.k2 <- par[4]; tau.rho <- par[5]
      
      # Constant term (normalization)
      const <- (1 - tau.rho^2) * (1 - tau.k1^2) * (1 - tau.k2^2)
      const <- (1 - tau.rho^2) * (1 - tau.k1^2) * (1 - tau.k2^2) / (4 * Pi^2)
      
      # Compute coefficients (from torus model definition)
      c0 <- (1 + tau.rho^2) * (1 + tau.k1^2) * (1 + tau.k2^2) - 8 * abs(tau.rho) * tau.k1 * tau.k2
      c1 <- 2 * (1 + tau.rho^2) * tau.k1 * (1 + tau.k2^2) - 4 * abs(tau.rho) * (1 + tau.k1^2) * tau.k2
      c2 <- 2 * (1 + tau.rho^2) * tau.k2 * (1 + tau.k1^2) - 4 * abs(tau.rho) * (1 + tau.k2^2) * tau.k1
      c3 <- -4 * (1 + tau.rho^2) * tau.k1 * tau.k2 + 2 * abs(tau.rho) * (1 + tau.k1^2) * (1 + tau.k2^2)
      c4 <- 2 * tau.rho * (1 - tau.k1^2) * (1 - tau.k2^2)
      
      # Evaluate torus density
      f <- const / (c0 - c1 * cos(theta1 - tau.mu1) - c2 * cos(theta2 - tau.mu2) -
                      c3 * cos(theta1 - tau.mu1) * cos(theta2 - tau.mu2) -
                      c4 * (sin(theta1 - tau.mu1) * sin(theta2 - tau.mu2)))
      return(f)
    },
    
    #' @description
    #' Computes the weighted negative log-likelihood for optimization.
    #' Used to re-estimate torus parameters given state responsibilities.
    #' 
    #' @param par Numeric vector of parameters (\eqn{\mu_1, \mu_2, \kappa_1, \kappa_2, \rho}).
    #' @param theta1,theta2 Numeric vectors of angle data.
    #' @param Pi Optional state probabilities.
    #' @param weights Numeric vector of posterior weights for each observation.
    #' @return Scalar negative weighted log-likelihood.
    weighted_loglikelihood = function(par, theta1, theta2, Pi, weights = 1) {
      # Transform bounded parameters for optimization stability
      par[3] <- (tanh(par[3]) + 1) / 2
      par[4] <- (tanh(par[4]) + 1) / 2
      par[5] <- tanh(par[5])
      
      # Evaluate density
      f <- self$density_function(theta1, theta2, Pi, par)
      val <- -sum(weights * log(f))
      return(val)
    },
    
    #' @description
    #' Updates the torus parameters by minimizing the weighted negative log-likelihood
    #' using numerical optimization (\code{optim}).
    #' 
    #' @param params List of model parameters (including data and number of states).
    #' @param Pi Optional vector of state probabilities.
    #' @param post.pi Matrix of posterior probabilities (if omitted, uses stored values).
    update = function(params, Pi = NULL, post.pi = NULL) {
      data <- params$data
      post.pi <- if (!is.null(post.pi)) post.pi else self$post.pi
      
      # Optimize parameters for each state
      for (k in 1:params$n_states) {
        init.par <- c(self$tau.mu1[k], self$tau.mu2[k],
                      self$tau.k1[k], self$tau.k2[k], self$tau.rho[k])
        
        opt <- optim(
          par = init.par,
          fn = self$weighted_loglikelihood,
          theta1 = data[, 1],
          theta2 = data[, 2],
          Pi = Pi,
          weights = post.pi[, k]
        )
        
        # Update parameters (with inverse transformations)
        self$tau.mu1[k] <- opt$par[1]
        self$tau.mu2[k] <- opt$par[2]
        self$tau.k1[k]  <- (tanh(opt$par[3]) + 1) / 2
        self$tau.k2[k]  <- (tanh(opt$par[4]) + 1) / 2
        self$tau.rho[k] <- tanh(opt$par[5])
      }
      
      # Recompute updated emission densities
      self$compute_density(params, Pi)
    }
  )
)
