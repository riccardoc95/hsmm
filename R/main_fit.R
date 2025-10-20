fit_hsmm <- function(
  #HMM
  data,
  n_states,
  covariates_omega = NULL,
  covariates_q = NULL,

  # HSMM
  semi = TRUE,
  max_dwell = NULL,
  
  # Model Type
  model_type = "torus",

  max_iter = 100,
  tol = 1e-5,
  verbose = TRUE,
  init = NULL,
  seed = NULL) {
  
  n_obs <- nrow(data)
  
  if(!is.null(max_dwell)){
    dwell_lengths <- rep(max_dwell, n_states)
    state_indices <- rep(1:n_states, dwell_lengths)
  } else {
    max_dwell <- 1
    dwell_lengths <- rep(max_dwell, n_states)
    state_indices <- rep(1:n_states, dwell_lengths)
  }

  params <- c(as.list(environment()))
  
  # duration model
  duration.model.factory <-DurationModelFactory$new()
  duration.model <- duration.model.factory$create(params)
  
  # transition model
  transition.model.factory <- TransitionModelFactory$new()
  transition.model <- transition.model.factory$create(params, 
                                                      duration.model$p.array)
  
  # emission model
  emission.model.factory <- EmissionModelFactory$new()
  emission.model <- emission.model.factory$create(params,
                                                  transition.model$Pi,
                                                  transition.model$post.pi)
  
  # EM-algorithm
  llk <- em.algorithm(params, 
                      emission.model, transition.model, duration.model)
  
  return(list(
    input_params = params,
    emission.model = emission.model,
    transition.model = transition.model,
    duration.model = duration.model
  ))
  
}
