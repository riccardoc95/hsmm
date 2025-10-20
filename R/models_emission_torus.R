library(R6)


TorusModel <- R6Class(
  "TorusModel",
  inherit = EmissionModel,
  
  public = list(
    tau.mu1 = NULL,
    tau.mu2 = NULL,
    tau.k1 = NULL,
    tau.k2 = NULL,
    tau.rho = NULL,
    density = NULL,
    
    initialize = function(params, Pi = NULL, post.pi = NULL) {
      super$initialize(params, Pi, post.pi)
      
      # Inizializza i parametri (tau.mu*, tau.k*, tau.rho)
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
      
      # Calcola la densità iniziale
      self$compute_density(params, Pi)
    },
    
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
    
    density_function = function(theta1, theta2, Pi, par = NULL) {
      tau.mu1 <- par[1]; tau.mu2 <- par[2]
      tau.k1 <- par[3]; tau.k2 <- par[4]; tau.rho <- par[5]
      
      const <- (1 - tau.rho^2) * (1 - tau.k1^2) * (1 - tau.k2^2)
      const <- (1 - tau.rho^2) * (1 - tau.k1^2) * (1 - tau.k2^2) / (4 * Pi^2)
      c0 <- (1 + tau.rho^2) * (1 + tau.k1^2) * (1 + tau.k2^2) - 8 * abs(tau.rho) * tau.k1 * tau.k2
      c1 <- 2 * (1 + tau.rho^2) * tau.k1 * (1 + tau.k2^2) - 4 * abs(tau.rho) * (1 + tau.k1^2) * tau.k2
      c2 <- 2 * (1 + tau.rho^2) * tau.k2 * (1 + tau.k1^2) - 4 * abs(tau.rho) * (1 + tau.k2^2) * tau.k1
      c3 <- -4 * (1 + tau.rho^2) * tau.k1 * tau.k2 + 2 * abs(tau.rho) * (1 + tau.k1^2) * (1 + tau.k2^2)
      c4 <- 2 * tau.rho * (1 - tau.k1^2) * (1 - tau.k2^2)
      
      f <- const / (c0 - c1 * cos(theta1 - tau.mu1) - c2 * cos(theta2 - tau.mu2) -
                      c3 * cos(theta1 - tau.mu1) * cos(theta2 - tau.mu2) -
                      c4 * (sin(theta1 - tau.mu1) * sin(theta2 - tau.mu2)))
      return(f)
    },
    
    weighted_loglikelihood = function(par, theta1, theta2, Pi, weights = 1) {
      par[3] <- (tanh(par[3])+1)/2
      par[4] <- (tanh(par[4])+1)/2
      par[5] <- tanh(par[5])
      
      f <- self$density_function(theta1, theta2, Pi, par)
      val <- -sum(weights * log(f))
      return(val)
    },
    
    update = function(params, Pi = NULL, post.pi = NULL) {
      data <- params$data
      post.pi <- if (!is.null(post.pi)) post.pi else self$post.pi
      
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
        
        self$tau.mu1[k] <- opt$par[1]
        self$tau.mu2[k] <- opt$par[2]
        self$tau.k1[k]  <- (tanh(opt$par[3]) + 1) / 2
        self$tau.k2[k]  <- (tanh(opt$par[4]) + 1) / 2
        self$tau.rho[k] <- tanh(opt$par[5])
      }
      
      # Ricalcola la densità aggiornata
      self$compute_density(params, Pi)
    }
  )
)