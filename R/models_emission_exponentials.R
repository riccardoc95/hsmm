#' Gaussian Emission Model for HMM/HSMM
#'
#' Models the emission distribution of each hidden state as a multivariate
#' Gaussian. Used within an HMM or HSMM framework.
#'
#' @field mu Matrix (n_states x n_dim) of mean vectors.
#' @field sigma 3D array (n_dim x n_dim x n_states) of covariance matrices.
#' @field density Matrix of computed emission densities for each observation and state.
#'
GaussianModel <- R6Class(
  "GaussianModel",
  inherit = EmissionModel,
  public = list(
    mu = NULL,
    sigma = NULL,
    density = NULL,
    
    #' @description
    #' Initializes Gaussian model parameters (\code{mu}, \code{sigma}) using weighted sample estimates.
    #' @param params List of model parameters including \code{data} and \code{n_states}.
    #' @param Pi Initial state probabilities.
    #' @param post.pi Posterior probabilities per state (from previous iteration).
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)
      
      # Initialize parameters based on data dimensionality
      n_dim <- ncol(params$data)
      mu_init <- matrix(0, params$n_states, n_dim)
      sigma_init <- array(0, dim = c(n_dim, n_dim, params$n_states))
      
      # Weighted estimation for each state's mean and covariance
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        sw <- sum(w)
        mu_init[k, ] <- sum(w * params$data) / sw
        centered <- sweep(params$data, 2, mu_init[k, ])
        centered <- as.matrix(centered)
        sigma_init[, , k] <- (t(centered) %*% (centered * w)) / sw
      }
      self$mu <- mu_init
      self$sigma <- sigma_init
      
      # Compute initial emission densities
      self$compute_density(params)
    },
    
    #' @description
    #' Computes Gaussian emission densities for all observations and states.
    compute_density = function(params, Pi) {
      dens <- array(dim = c(params$n_obs, params$n_states))
      for (k in 1:params$n_states) {
        dens[, k] <- mvtnorm::dmvnorm(params$data, mean = self$mu[k, ], sigma = self$sigma[, , k])
      }
      self$density <- dens
    },
    
    #' @description
    #' Updates the Gaussian model parameters (\code{mu}, \code{sigma}) using new posterior weights.
    update = function(params, Pi, post.pi) {
      data <- params$data
      n_dim <- ncol(data)
      mu_new <- matrix(0, params$n_states, n_dim)
      sigma_new <- array(0, dim = c(n_dim, n_dim, params$n_states))
      
      # Re-estimate weighted mean and covariance
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        sw <- sum(w)
        mu_new[k, ] <- sum(w * params$data) / sw
        centered <- sweep(data, 2, mu_new[k, ])
        centered <- as.matrix(centered)
        sigma_new[, , k] <- (t(centered) %*% (centered * w)) / sw
      }
      self$mu <- mu_new
      self$sigma <- sigma_new
      
      self$compute_density(params)
    }
  )
)


#' Poisson Emission Model
#'
#' Models emission probabilities using a Poisson distribution for count data.
#'
#' @field lambda Numeric vector of Poisson rate parameters (\eqn{\lambda_k}).
#' @field density Matrix of emission probabilities for each observation and state.
#'
PoissonModel <- R6Class(
  "PoissonModel",
  inherit = EmissionModel,
  public = list(
    lambda = NULL,
    n_states = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)
      
      # Weighted initialization of lambda for each state
      lambda_init <- numeric(params$n_states)
      data <- params$data
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        lambda_init[k] <- sum(w * data) / sum(w)
      }
      self$lambda <- lambda_init
      
      self$compute_density(params)
      invisible(self)
    },
    
    compute_density = function(params, Pi) {
      data <- params$data
      n_obs <- nrow(data)
      dens <- array(dim = c(n_obs, params$n_states))
      for (k in 1:params$n_states) {
        dens[, k] <- dpois(data, lambda = self$lambda[k])
      }
      self$density <- dens
      invisible(self)
    },
    
    update = function(params, Pi, post.pi) {
      data <- params$data
      lambda_new <- numeric(params$n_states)
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        lambda_new[k] <- sum(w * data) / sum(w)
      }
      self$lambda <- lambda_new
      self$compute_density(params)
      invisible(self)
    }
  )
)


#' Exponential Emission Model
#'
#' Models continuous non-negative data using an exponential distribution.
#'
#' @field lambda Numeric vector of rate parameters.
#' @field density Matrix of emission probabilities for each observation and state.
#'
ExponentialModel <- R6Class(
  "ExponentialModel",
  inherit = EmissionModel,
  public = list(
    lambda = NULL,
    n_states = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)
      
      # Weighted initialization of lambda for each state
      lambda_init <- numeric(params$n_states)
      data <- params$data
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        lambda_init[k] <- sum(w) / sum(w * data)
      }
      self$lambda <- lambda_init
      self$compute_density(params)
      invisible(self)
    },
    
    compute_density = function(params, Pi) {
      data <- params$data
      n_obs <- nrow(data)
      dens <- array(dim = c(n_obs, params$n_states))
      for (k in 1:params$n_states) {
        dens[, k] <- dexp(data, rate = self$lambda[k])
      }
      self$density <- dens
      invisible(self)
    },
    
    update = function(params, Pi, post.pi) {
      data <- params$data
      lambda_new <- numeric(params$n_states)
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        lambda_new[k] <- sum(w) / sum(w * data)
      }
      self$lambda <- lambda_new
      self$compute_density(params)
      invisible(self)
    }
  )
)


#' Gamma Emission Model
#'
#' Models emission probabilities using a Gamma distribution.
#'
#' @field shape Numeric vector of shape parameters.
#' @field scale Numeric vector of scale parameters.
#' @field density Matrix of emission probabilities for each observation and state.
#'
GammaModel <- R6Class(
  "GammaModel",
  inherit = EmissionModel,
  public = list(
    shape = NULL,
    scale = NULL,
    n_states = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)
      
      # Weighted initialization of shape and scale
      shape_init <- numeric(params$n_states)
      scale_init <- numeric(params$n_states)
      data <- params$data
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        m <- sum(w * data) / sum(w)
        v <- sum(w * (data - m)^2) / sum(w)
        shape_init[k] <- m^2 / v
        scale_init[k] <- v / m
      }
      self$shape <- shape_init
      self$scale <- scale_init
      self$compute_density(params)
      invisible(self)
    },
    
    compute_density = function(params, Pi) {
      data <- params$data
      n_obs <- nrow(data)
      dens <- array(dim = c(n_obs, params$n_states))
      for (k in 1:params$n_states) {
        dens[, k] <- dgamma(data, shape = self$shape[k], scale = self$scale[k])
      }
      self$density <- dens
      invisible(self)
    },
    
    update = function(params, Pi, post.pi) {
      data <- params$data
      shape_new <- numeric(params$n_states)
      scale_new <- numeric(params$n_states)
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        m <- sum(w * data) / sum(w)
        v <- sum(w * (data - m)^2) / sum(w)
        shape_new[k] <- m^2 / v
        scale_new[k] <- v / m
      }
      self$shape <- shape_new
      self$scale <- scale_new
      self$compute_density(params)
      invisible(self)
    }
  )
)


#' Beta Emission Model
#'
#' Models emission probabilities for variables bounded in [0, 1] using the Beta distribution.
#'
#' @field alpha Shape parameter α for each state.
#' @field beta Shape parameter β for each state.
#' @field density Matrix of emission probabilities for each observation and state.
#'
BetaModel <- R6Class(
  "BetaModel",
  inherit = EmissionModel,
  public = list(
    alpha = NULL,
    beta = NULL,
    n_states = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)
      
      # Weighted initialization of alpha and beta parameters
      alpha_init <- numeric(params$n_states)
      beta_init <- numeric(params$n_states)
      data <- params$data
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        m <- sum(w * data) / sum(w)
        v <- sum(w * (data - m)^2) / sum(w)
        tmp <- m * (1 - m) / v - 1
        alpha_init[k] <- m * tmp
        beta_init[k] <- (1 - m) * tmp
      }
      self$alpha <- alpha_init
      self$beta <- beta_init
      self$compute_density(params)
      invisible(self)
    },
    
    compute_density = function(params, Pi) {
      data <- params$data
      n_obs <- nrow(data)
      dens <- array(dim = c(n_obs, params$n_states))
      for (k in 1:params$n_states) {
        dens[, k] <- dbeta(data, shape1 = self$alpha[k], shape2 = self$beta[k])
      }
      self$density <- dens
      invisible(self)
    },
    
    update = function(params, Pi, post.pi) {
      data <- params$data
      alpha_new <- numeric(params$n_states)
      beta_new <- numeric(params$n_states)
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        m <- sum(w * data) / sum(w)
        v <- sum(w * (data - m)^2) / sum(w)
        tmp <- m * (1 - m) / v - 1
        alpha_new[k] <- m * tmp
        beta_new[k] <- (1 - m) * tmp
      }
      self$alpha <- alpha_new
      self$beta <- beta_new
      self$compute_density(params)
      invisible(self)
    }
  )
)


#' Student-t Emission Model
#'
#' Models emission probabilities using a Student-t distribution, allowing for
#' heavier tails compared to the Gaussian model.
#'
#' @field mu Mean for each state.
#' @field sigma Scale (standard deviation) for each state.
#' @field df Degrees of freedom.
#' @field density Matrix of emission probabilities for each observation and state.
#'
StudentTModel <- R6Class(
  "StudentTModel",
  inherit = EmissionModel,
  public = list(
    mu = NULL,
    sigma = NULL,
    df = NULL,
    n_states = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)
      
      # Weighted initialization for mu and sigma
      mu_init <- numeric(params$n_states)
      sigma_init <- numeric(params$n_states)
      df_init <- rep(10, params$n_states)  # Default degrees of freedom
      
      data <- params$data
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        mu_init[k] <- sum(w * data) / sum(w)
        sigma_init[k] <- sqrt(sum(w * (data - mu_init[k])^2) / sum(w))
      }
      self$mu <- mu_init
      self$sigma <- sigma_init
      self$df <- df_init
      self$compute_density(params)
      invisible(self)
    },
    
    compute_density = function(params, Pi) {
      data <- params$data
      n_obs <- nrow(data)
      dens <- array(dim = c(n_obs, params$n_states))
      for (k in 1:params$n_states)
        dens[, k] <- dt((data - self$mu[k]) / self$sigma[k], df = self$df[k]) / self$sigma[k]
      self$density <- dens
      invisible(self)
    },
    
    update = function(params, Pi, post.pi) {
      data <- params$data
      mu_new <- numeric(params$n_states)
      sigma_new <- numeric(params$n_states)
      df_new <- self$df  # Keep df fixed by default
      
      # Weighted updates
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        mu_new[k] <- sum(w * data) / sum(w)
        sigma_new[k] <- sqrt(sum(w * (data - mu_new[k])^2) / sum(w))
      }
      self$mu <- mu_new
      self$sigma <- sigma_new
      self$df <- df_new
      self$compute_density(params)
      invisible(self)
    }
  )
)
