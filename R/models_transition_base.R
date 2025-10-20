library(R6)


TransitionModel <- R6Class(
  "TransitionModel",
  public = list(
    gamma = NULL,
    omega = NULL,
    Pi = NULL,
    post.pi = NULL,
    initialize = function(params) {
      # pi
      self$Pi <- runif(params$n_states)
      self$Pi <- self$Pi / sum(self$Pi)
      
      # post.pi
      self$post.pi <- LaplacesDemon::rdirichlet(params$n_obs, alpha = rep(1 /
                                                                            params$n_states, params$n_states))
      
    },
    compute_gamma = function(params, p.array) {
      invisible(self)
    },
    compute_omega = function(params, post.bi.pi.aug) {
      invisible(self)
    },
    update = function(params, post.pi.aug){
      self$post.pi <- t(sapply(1:params$n_obs, 
                                           function(i){sapply(split(post.pi.aug[i,], 
                                                                    params$state_indices), 
                                                              sum)}
                                           )
                                    )
      self$Pi <- self$post.pi[1,]
    }
  )
)
    