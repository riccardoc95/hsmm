GaussianModel <- R6Class(
  "GaussianModel",
  inherit = EmissionModel,
  public = list(
    mu = NULL,
    sigma = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)

      # Inizializzazione parametri
      n_dim <- ncol(params$data)
      mu_init <- matrix(0, params$n_states, n_dim)
      sigma_init <- array(0, dim = c(n_dim, n_dim, params$n_states))
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
      
      # Calcolo densità iniziale
      self$compute_density(params)
    },
    
    compute_density = function(params, Pi) {
      dens <- array(dim = c(params$n_obs, params$n_states))
      for (k in 1:params$n_states) {
        dens[, k] <- mvtnorm::dmvnorm(params$data, mean = self$mu[k, ], sigma = self$sigma[, , k])
      }
      self$density <- dens
    },
    
    update = function(params, Pi, post.pi) {
      data <- params$data
      n_dim <- ncol(data)
      mu_new <- matrix(0, params$n_states, n_dim)
      sigma_new <- array(0, dim = c(n_dim, n_dim, params$n_states))
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



PoissonModel <- R6Class(
  "PoissonModel",
  inherit = EmissionModel,
  public = list(
    lambda = NULL,
    n_states = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)

      # Inizializzazione lambda
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


ExponentialModel <- R6Class(
  "ExponentialModel",
  inherit = EmissionModel,
  public = list(
    lambda = NULL,
    n_states = NULL,
    density = NULL,
    
    initialize = function(params, Pi, post.pi) {
      super$initialize(params, Pi, post.pi)
      
      
      # Inizializzazione lambda
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
      
      
      # Inizializzazione shape e scale
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
      
      
      # Inizializzazione alpha e beta
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
      
      
      # Inizializzazione parametri mu, sigma, df
      mu_init <- numeric(params$n_states)
      sigma_init <- numeric(params$n_states)
      # df iniziale fisso, es 10
      df_init <- rep(10, params$n_states)
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
      df_new <- self$df  # manteniamo df fisso o potresti definire altra logica
      
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
