em.step = function(params, emission.model, transition.model, duration.model){
  # E-step
  density <- emission.model$density
  
  fit <- matrix(0, nrow=nrow(density), ncol=params$n_states*params$max_dwell)
  for(k in 1:params$n_states){
    pos <- ((k-1)*params$max_dwell + 1):(k*params$max_dwell)
    fit[, pos] <- matrix(rep(density[, k], params$max_dwell), 
                         nrow=nrow(density), ncol=params$max_dwell)
  }
  
  # Backward forward algorithm
  bw <- backward.forward(fit, transition.model$gamma, transition.model$Pi, 
                         params$n_states, params$max_dwell)
  #bw <- lbackward_forward(fit, transition.model$gamma, transition.model$Pi, 
  #                        params$n_states, params$max_dwell)

  
  # M-step
  
  # Update post.pi/Pi
  transition.model$update(params, bw$post.pi.aug )
  
  # First: maximization of the Omegas
  transition.model$compute_omega(params, bw$post.bi.pi.aug)
  
  # Second: parametric estimation of the dwell-times
  duration.model$compute_beta(params, bw$post.bi.pi.aug)
  
  
  # Update p.array
  duration.model$compute_p.array(params)
  
  # Update Gamma
  transition.model$compute_gamma(params, duration.model$p.array)
  
  # Update the density parameters
  emission.model$update(params, 
                        transition.model$Pi, 
                        transition.model$post.pi)
  
  return(bw$l)
}

em.algorithm = function(params, 
                        emission.model, 
                        transition.model, 
                        duration.model){
  step <- 1
  llk <- numeric(params$max_iter)
  converged <- FALSE
  dif <- NA
  
  while ( !converged && step <= params$max_iter) {
    llk[step] <- em.step(params, 
                         emission.model, transition.model, duration.model)
    
    if(step > 10){
      dif <- abs(llk[step]-old.llk)
      if (!is.na(dif) && dif < params$tol){
        converged <- TRUE
      }
    }
    
    old.llk <- llk[step]
    step <- step + 1
    
    if(params$verbose){
      cat("iteration ", step-1,"; loglik = ", llk[step-1], "\n")
    } 
  }
  return(llk)
}