library(R6)

#############################
# TODO: NOT WORK!!!!!!!!!!  #
#############################

VarGaussianModel <- R6::R6Class(
  classname = "VarGaussianModel",
  inherit = EmissionModel,
  public = list(
    auto = NULL,       
    Sigma = NULL,      
    ar = NULL,        
    y = NULL,         
    lambda = 0,        
    fitted = FALSE,
    
    initialize = function(params, Pi, post.pi) {
      # dati
      self$y <- params$data
      
      self$ar <- if (is.null(params$ar)) 1 else params$ar
      self$lambda <- if (is.null(params$lambda)) 0 else params$lambda
      
      K <- params$n_states
      M <- params$max_dwell
      n <- nrow(self$data)
      p <- ncol(self$data)
      
      if (self$ar > 0) {
        ylag <- self$y[1:(n - self$ar), ]
        for (lag in 2:self$ar) {
          ylag <- cbind(ylag, self$y[lag:(n - self$ar + lag - 1), ])
        }
        colnames(ylag) <- paste0("lag", 1:ncol(ylag))
        yout <- self$y[(1 + self$ar):n, ]
      } else {
        stop("VAR Gaussian requires ar >= 1")
      }
      
      post.pi <- as.matrix(post.pi)
      J <- ncol(ylag) + 1
      self$auto <- vector("list", K)
      self$Sigma <- vector("list", K)
      
      for (k in 1:K) {
        self$auto[[k]] <- matrix(NA, nrow = J, ncol = p)
        lambdak <- self$lambda * sqrt(sum(post.pi[, k]))
        
        for (j in 1:p) {
          fit_glmnet <- glmnet::glmnet(
            x = ylag,
            y = yout[, j],
            family = "gaussian",
            weights = post.pi[, k],
            lambda = lambdak,
            standardize = TRUE
          )
          self$auto[[k]][, j] <- c(fit_glmnet$a0, fit_glmnet$beta[, 1])
        }
        
        resid <- yout - cbind(1, ylag) %*% self$auto[[k]]
        self$Sigma[[k]] <- t((resid * post.pi[, k]) / sum(post.pi[, k])) %*% resid
      }
      
      self$fitted <- TRUE
      invisible(self)
    },
    
    compute_density = function(params, p.array) {
      if (!self$fitted)
        stop("Model not initialized or fitted yet.")
      
      y <- self$y
      n <- nrow(y)
      p <- ncol(y)
      ar <- self$ar
      K <- params$n_states
      M <- params$max_dwell
      ld <- rep(M, K)
      d <- rep(1:K, ld)
      yout <- y[(1 + ar):n, ]
      ylag <- y[1:(n - ar), , drop = FALSE]
      for (lag in 2:ar) {
        ylag <- cbind(ylag, y[lag:(n - ar + lag - 1), , drop = FALSE])
      }
      
      intxwlag <- cbind(1, ylag)
      fit <- matrix(0, nrow = n - ar, ncol = sum(ld))
      for (k in 1:K) {
        pos <- which(d == k)
        mus <- intxwlag %*% self$auto[[k]]
        for (i in 1:(n - ar)) {
          fit[i, pos] <- mvtnorm::dmvnorm(yout[i, ], mean = mus[i, ], sigma = self$Sigma[[k]])
        }
      }
      
      self$density <- fit
      invisible(self)
    },
    

    update = function(params, Pi, post.pi) {
      if (is.null(self$y))
        stop("Data missing in VAR Gaussian emission model.")
      self$initialize(params, Pi, post.pi)
      invisible(self)
    }
    
  )
)
