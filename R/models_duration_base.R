library(R6)

DurationModel <- R6Class(
  "DurationModel",
  public = list(
    #p.array
    p.array = NULL,
    # beta
    beta.intercepts = NULL,
    # c'è sempre
    beta.time_effects = NULL,
    # semi = TRUE -> anche questa!!
    beta.covariate = NULL,
    # q_variance = NON VUOTO -> anche questa!!
    
    d.state = NULL,
    # (d_matrix) # c'è sempre
    d.time = NULL,
    # (d_matrix) # semi = TRUE -> anche questa!!
    d.covariates = NULL,
    # (d_matrix) # q_variance = NON VUOTO -> anche questa!!
    
    link = NULL,
    initialize = function(params) {
      self$beta.intercepts <- runif(params$n_states, -5, 0)
      self$d.state <- rep(1:params$n_states, each=params$max_dwell)
      if (params$semi){
        self$beta.time_effects <- runif(params$n_states, 0.01, 0.4)
        self$d.time <- rep(1:params$max_dwell, times=params$n_states) + 0.5
      } else {
        self$beta.time_effects <- rep(0, length(self$beta.intercepts))
        self$d.time <- rep(0, length(self$d.state))
      }
      
      self$link <- "cloglog"
    },
    
    gomp = function(i, alpha, betai){
      exp(alpha + betai*i)
    },
    
    compute_p.array = function(params) {
      invisible(self)
    },
    
    compute_beta = function(params, post.bi.pi.aug) {
      invisible(self)
    }
  )
)