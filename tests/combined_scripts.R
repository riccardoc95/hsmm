# ---- Start of algorithms_bf.R ----


backward.forward <- function(fit, Gamma, Pi, K, M){
  
  ld <- rep(M, K)
  n <- nrow(fit)
  
  psi <- alpha <- Beta <- array(dim = c(n, sum(ld)))
  C	<- array(dim = n)
  
  post.pi.aug <- array(dim = c(n, sum(ld)))
  post.bi.pi.aug <- array(dim = c(n-1, sum(ld), sum(ld)))
  
  psi[1,] <- Pi
  psi_per_fit <- psi[1,]*fit[1,]
  C[1] <- sum(psi_per_fit)
  alpha[1,] <- (psi_per_fit)/C[1] 
  
  for(i in 2:n){
    psi[i, ] <- alpha[i-1,]%*%Gamma[i-1,,]
    psi_per_fit <- psi[i,]*fit[i,]
    C[i] <- sum(psi_per_fit)
    alpha[i,] <- (psi_per_fit)/C[i] 
    
  }
  
  l <- sum(log(C))
  
  Beta[n,] <- 1/C[n]
  
  for(i in (n-1):1){
    Beta[i, ] <- c(Gamma[i,,]%*%(fit[i+1,]*Beta[i+1,]))/C[i]
  }
  
  alphaBeta <- alpha*Beta
  post.pi.aug = alphaBeta/rowSums(alphaBeta)
  
  for(h in 1:sum(ld)){
    for(k in 1:sum(ld)){
      post.bi.pi.aug[,h,k] <- alpha[1:(n-1),h]*Gamma[,h,k]*fit[2:n,k]*Beta[2:n,k] 
    }
  }
  return(list(l = l, post.pi.aug = post.pi.aug, post.bi.pi.aug = post.bi.pi.aug))
}

# ---- End of algorithms_bf.R ----

# ---- Start of algorithms_em.R ----

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

# ---- End of algorithms_em.R ----

# ---- Start of algorithms_optim.R ----



# ---- End of algorithms_optim.R ----

# ---- Start of main_fit.R ----

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

# ---- End of main_fit.R ----

# ---- Start of models_duration_base.R ----

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

# ---- End of models_duration_base.R ----

# ---- Start of models_duration_geometric.R ----

library(R6)

DurationModel.Geometric <- R6Class(
  "DurationModel.Geometric",
  inherit = DurationModel,
  public = list(
    initialize = function(params) {
      super$initialize(params)
      
      self$beta.covariate <- NULL
      self$d.covariates <- NULL
      self$link <- "cloglog"
      self$compute_p.array(params)
    },
    
    p.k = function(i, alpha, betai){
      1 - exp(-(self$gomp(i + 0.5, alpha, betai)))
    }, 
    
    compute_p.array = function(params) {
      self$p.array <- array(0, 
                            dim = c(params$n_obs-1, sum(params$dwell_lengths))) 
      
      for(t in 1:(params$n_obs-1)){ 
        for(k in 1:params$n_states){
          pos <- which(params$state_indices == k)
          self$p.array[t, pos] <- self$p.k(1:params$dwell_lengths[k], 
                                      alpha = self$beta.intercepts[k], 
                                      betai = self$beta.time_effects[k])
        }
      }
    },
    
    compute_beta = function(params, post.bi.pi.aug) {
      n <- params$n_obs - 1
      d <- rep(1:params$n_states, params$dwell_lengths)
      
      cases_fin <- noncases_fin <- numeric()
      for (t in 1:n) {
        for (h in 1:params$n_states) {
          posh <- which(d == h)
          postbipi_1 <- post.bi.pi.aug[t, posh, -posh, drop = FALSE]
          postbipi_2 <- post.bi.pi.aug[t, posh, posh, drop = FALSE]
          cases_fin <- c(cases_fin, rowSums(postbipi_1))
          noncases_fin <- c(noncases_fin, rowSums(postbipi_2))
        }
      }
      
      # Costruisci d_matrix con solo state e time, niente covariate
      d_matrix <- data.frame(state = self$d.state, time = self$d.time)
      
      outno0 <- (cases_fin + noncases_fin) != 0
      df_for_glm <- data.frame(cases = cases_fin[outno0], noncases = noncases_fin[outno0], d_matrix[outno0, ])
      
      formula_dwell <- as.formula("cbind(cases, noncases) ~ factor(state) * time")
      
      fit_beta <- suppressWarnings(glm(formula_dwell, family = binomial(link = self$link), data = df_for_glm))
      
      coefs <- coef(fit_beta)
      
      # Parsing coefficienti solo per stato e tempo
      intercept_ids <- grep("^factor\\(state\\)[0-9]+$", names(coefs))
      time_id <- grep("time", names(coefs))
      
      beta_intercepts <- numeric(params$n_states)
      beta_intercepts[1] <- coefs["(Intercept)"]
      if(length(intercept_ids) > 0) {
        for(i in 2:params$n_states) {
          nm <- paste0("factor(state)", i)
          if(nm %in% names(coefs)) {
            beta_intercepts[i] <- beta_intercepts[1] + coefs[nm]
          } else {
            beta_intercepts[i] <- beta_intercepts[1]
          }
        }
      } else {
        beta_intercepts[] <- coefs["(Intercept)"]
      }
      
      beta_time_effects <- numeric(params$n_states)
      beta_time_effects[] <- coefs[time_id][1]
      
      # Aggiorna beta
      self$beta.intercepts <- beta_intercepts
      if (params$semi){
        self$beta.time_effects <- beta_time_effects
      } else {
        self$beta.time_effects <- rep(0, length(self$beta.intercepts))
      }
      self$beta.covariate <- NULL  
      
    }
  )
)

# ---- End of models_duration_geometric.R ----

# ---- Start of models_duration_gompertz.R ----

library(R6)

DurationModel.Gompertz <- R6Class(
  "DurationModel.Gompertz",
  inherit = DurationModel,
  public = list(
    initialize = function(params) {
      super$initialize(params)
      
      self$beta.covariate <-  array(
        runif(params$n_states * ncol(params$covariates_q), -0.5, 0.5),
        dim = c(params$n_states, ncol(params$covariates_q)))
      
      self$d.covariates <- matrix(
        NA,
        nrow = sum(params$dwell_lengths) * (params$n_obs - 1),
        ncol = ncol(params$covariates_q)
      )
      for (j in 1:ncol(params$covariates_q)) {
        self$d.covariates[, j] <- rep(params$covariates_q[1:(params$n_obs - 1), j], 
                                      each = sum(params$dwell_lengths))
      }
      
      self$link <- "cloglog"
      self$compute_p.array(params)
    },
    
    p.k = function(i, alpha, betai, beta_x, x_t){
      1 - exp(
        -(
          self$gomp(i + 0.5, alpha, betai)*exp(as.numeric(
            # TODO: sistemare!!
            drop(as.matrix(x_t)%*%as.matrix(beta_x)))
            )
          )
        )
    }, 
    
    compute_p.array = function(params) {
      self$p.array <- array(0, 
                            dim = c(params$n_obs-1, sum(params$dwell_lengths))) 
      
      for(t in 1:(params$n_obs-1)){ 
        for(k in 1:params$n_states){
          pos <- which(params$state_indices == k)
          self$p.array[t, pos] <- self$p.k(1:params$dwell_lengths[k], 
                                           alpha = self$beta.intercepts[k], 
                                           betai = self$beta.time_effects[k], 
                                           beta_x = self$beta.covariate[k,], 
                                           x_t = params$covariates_q[t, ])
        }
      }
    },
    
    compute_beta = function(params, post.bi.pi.aug) {
      n <- params$n_obs - 1
      d <- rep(1:params$n_states, params$dwell_lengths)
      q <- ncol(params$covariates_q)
      
      cases_fin <- noncases_fin <- numeric()
      for (t in 1:n) {
        for (h in 1:params$n_states) {
          posh <- which(d == h)
          postbipi_1 <- post.bi.pi.aug[t, posh, -posh, drop = FALSE]
          postbipi_2 <- post.bi.pi.aug[t, posh, posh, drop = FALSE]
          cases_fin <- c(cases_fin, rowSums(postbipi_1))
          noncases_fin <- c(noncases_fin, rowSums(postbipi_2))
        }
      }
      
      # Costruisci d_matrix con colonne: State, Time, Covariate1,...,Covariate q
      d_matrix <-  data.frame(state = self$d.state, 
                              time = self$d.time,
                              setNames(as.data.frame(self$d.covariates), 
                                       paste0("V", 1:ncol(self$d.covariates))))
  
      outno0 <- (cases_fin + noncases_fin) != 0
      df_for_glm <- data.frame(cases = cases_fin[outno0], noncases = noncases_fin[outno0], d_matrix[outno0, ])
      formula_str <- paste0("cbind(cases, noncases) ~ factor(state)*time + ",
                            paste0(paste0("factor(state)*V", 1:q), collapse = " + "))
      formula_dwell <- as.formula(formula_str)
      
      fit_beta <- suppressWarnings(glm(formula_dwell, family = binomial(link = self$link), data = df_for_glm))
      
      coefs <- coef(fit_beta)
      
      
      # Intercetta e coefficienti per stato (es. factor(state)2, factor(state)3, ...)
      intercept_ids <- grep("^factor\\(state\\)[0-9]+$", names(coefs))
      # Coefficiente per tempo (es. time)
      time_id <- grep("time", names(coefs))
      # Coefficienti per covariate per ogni stato (es. factor(state)2:V1, factor(state)3:V2, ...)
      covar_ids <- grep("V", names(coefs))
      
      # Aggiorna beta.intercepts
      beta_intercepts <- numeric(params$n_states)
      beta_intercepts[1] <- coefs["(Intercept)"]
      if(length(intercept_ids) > 0) {
        for(i in 2:params$n_states) {
          nm <- paste0("factor(state)", i)
          if(nm %in% names(coefs)) {
            beta_intercepts[i] <- beta_intercepts[1] + coefs[nm]
          } else {
            beta_intercepts[i] <- beta_intercepts[1]
          }
        }
      } else {
        beta_intercepts[] <- coefs["(Intercept)"]
      }
      
      # Aggiorna beta.time_effects (assumendo effetto tempo comune o per stato)
      beta_time_effects <- numeric(params$n_states)
      beta_time_effects[] <- coefs[time_id][1]  # se unico coefficiente per tempo
      
      # Aggiorna beta.covariate matrice n_states x q
      q <- length(covar_ids) / params$n_states
      beta_covar <- matrix(0, nrow = params$n_states, ncol = q)
      for (j in 1:q) {
        base_name <- paste0("V", j)
        beta_covar[1, j] <- coefs[grep(paste0("^", base_name, "$"), names(coefs))]
        for (k in 2:params$n_states) {
          nm <- paste0("factor(state)", k, ":", base_name)
          if(nm %in% names(coefs)) {
            beta_covar[k, j] <- beta_covar[1, j] + coefs[nm]
          } else {
            beta_covar[k, j] <- beta_covar[1, j]
          }
        }
      }
      
      # Assegna
      self$beta.intercepts <- beta_intercepts
      if (params$semi){
        self$beta.time_effects <- beta_time_effects
      } else {
        self$beta.time_effects <- rep(0, length(self$beta.intercepts))
      }
      self$beta.covariate <- beta_covar
      
    }
  )
)

# ---- End of models_duration_gompertz.R ----

# ---- Start of models_emission_base.R ----

library(R6)


EmissionModel <- R6Class(
  "EmissionModel",
  public = list(
    density = NULL,
    initialize = function(params, Pi, post.pi) {
      invisible(self)
    },
    compute_density = function(params, p.array) {
      invisible(self)
    },
    update = function(params, Pi, post.pi) {
      invisible(self)
    }
  )
)

# ---- End of models_emission_base.R ----

# ---- Start of models_emission_exponentials.R ----

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
        mu_init[k, ] <- colSums(w * params$data) / sw
        centered <- sweep(params$data, 2, mu_init[k, ])
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
      post.pi <- params$post.pi
      n_dim <- ncol(data)
      mu_new <- matrix(0, params$n_states, n_dim)
      sigma_new <- array(0, dim = c(n_dim, n_dim, params$n_states))
      for (k in 1:params$n_states) {
        w <- post.pi[, k]
        sw <- sum(w)
        mu_new[k, ] <- colSums(w * data) / sw
        centered <- sweep(data, 2, mu_new[k, ])
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
      post.pi <- params$post.pi
      params$n_states <- ncol(post.pi)
      
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
      post.pi <- params$post.pi
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
      post.pi <- params$post.pi
      params$n_states <- ncol(post.pi)
      
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
      post.pi <- params$post.pi
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
      post.pi <- params$post.pi
      params$n_states <- ncol(post.pi)
      
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
      post.pi <- params$post.pi
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
      post.pi <- params$post.pi
      params$n_states <- ncol(post.pi)
      
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
      post.pi <- params$post.pi
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
      post.pi <- params$post.pi
      params$n_states <- ncol(post.pi)
      
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
      post.pi <- params$post.pi
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

# ---- End of models_emission_exponentials.R ----

# ---- Start of models_emission_torus.R ----

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
      safe_log <- function(x, eps = 1e-10) {
        log(pmax(x, eps))
      }
      
      f <- self$density_function(theta1, theta2, Pi, par)
      val <- -sum(weights * safe_log(f))
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

# ---- End of models_emission_torus.R ----

# ---- Start of models_factory.R ----

DurationModelFactory <- R6Class("DurationModelFactory", public = list(
  create = function(params) {
    if (is.null(params$covariates_q)){
      type <- "geometric"
    } else{
      type <- "gompertz"
    }
    
    switch(
      type,
      "geometric" = DurationModel.Geometric$new(params),
      "gompertz" = DurationModel.Gompertz$new(params),
      stop("Unknown transition model type: ", type)
    )
  }
))

TransitionModelFactory <- R6Class("TransitionModelFactory", public = list(
  create = function(params, p.array) {
    if (is.null(params$covariates_omega)){
      type <- "homogeneous"
    } else{
      type <- "timevaring"
    }
    switch(
      type,
      "homogeneous" = TransitionModel.Homogeneous$new(params, p.array),
      "timevaring" = TransitionModel.TimeVarying$new(params, p.array),
      stop("Unknown transition model type: ", type)
    )
  }
))

EmissionModelFactory <- R6Class("EmissionModelFactory", public = list(
  create = function(params, Pi, post.pi) {
    type <- params$model_type
    switch(
      type,
      "torus" = TorusModel$new(params, Pi, post.pi),
      "gaussian" = GaussianModel$new(params, Pi, post.pi),
      "poisson" = PoissonModel$new(params, Pi, post.pi),
      "exponential" = ExponentialModel$new(params, Pi, post.pi),
      "gamma" = GammaModel$new(params, Pi, post.pi),
      "beta" = BetaModel$new(params, Pi, post.pi),
      "student" = StudentTModel$new(params, Pi, post.pi),
      stop("Unknown emission model type: ", type)
    )
  }
))

# ---- End of models_factory.R ----

# ---- Start of models_transition_base.R ----

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
    

# ---- End of models_transition_base.R ----

# ---- Start of models_transition_homogeneous.R ----

library(R6)

TransitionModel.Homogeneous <- R6Class(
  "TransitionModel.Homogeneous",
  inherit = TransitionModel,
  public = list(
    initialize = function(params, p.array) {
      super$initialize(params)
      # omega
      if (params$n_states == 2) {
        self$omega <- matrix(c(0, 1, 1, 0), 2, 2)
      } else {
        self$omega <- matrix(runif(params$n_states * params$n_states),
                             params$n_states,
                             params$n_states)
        diag(self$omega) <- 0
        for (k in 1:params$n_states) {
          self$omega[k, -k] <- self$omega[k, -k] / sum(self$omega[k, -k])
        }
      }
      
      # gamma
      self$gamma <- self$compute_gamma(params, p.array)
      
    },
    compute_gamma = function(params, p.array) {
      K <- params$n_states
      M <- params$max_dwell
      Omega <- self$omega
      
      # function(Omega, p.array, K, M) {
      # Old code
      n <- nrow(p.array)
      ld <- rep(M, K)
      d <- rep(1:K, ld)
      p <- p.array
      Gamma <- array(0, dim = c(n, sum(ld), sum(ld)))
      for (t in 1:n) {
        for (k in 1:K) {
          pos <- which(d == k)
          if (length(pos) >= 2) {
            if (length(pos) == 2) {
              Gamma[t, min(pos):(max(pos) - 1), (min(pos) + 1):max(pos)] <- 1 - p[t, pos][-length(p[t, pos])]
            } else{
              diag(Gamma[t, min(pos):(max(pos) - 1), (min(pos) + 1):max(pos)]) <- 1 -
                p[t, pos][-length(p[t, pos])]
            }
          }
          Gamma[t, max(pos), max(pos)] <- 1 - p[t, pos][length(p[t, pos])]
          
          for (i in (1:K)[-k]) {
            Gamma[t, min(pos):max(pos), c(1, which(diff(d) == 1) + 1)[i]] <- p[t, pos] *
              Omega[k, i]
          }
        }
        
      }
      
      self$gamma <- Gamma
    },
    compute_omega = function(params, post.bi.pi.aug) {
      if (params$n_states > 2) {
        self$omega <- matrix(0, params$n_states, params$n_states)
        for (h in 1:params$n_states) {
          for (k in setdiff(1:params$n_states, h)) {
            self$omega[h, k] <- sum(post.bi.pi.aug[, h, k])
          }
          total <- sum(self$omega[h, -h])
          if (total > 0) {
            self$omega[h, -h] <- self$omega[h, -h] / total
          }
          self$omega[h, h] <- 0
        }
      }
    }
  )
)

# ---- End of models_transition_homogeneous.R ----

# ---- Start of models_transition_time_varying.R ----

library(R6)

TransitionModel.TimeVarying <- R6Class(
  "TransitionModel.TimeVarying",
  inherit = TransitionModel,
  public = list(
    initialize = function(params, p.array) {
      super$initialize(params)
      # omega
      if (params$n_states == 2){
        self$omega <- array(0, dim = c(params$n_obs - 1, params$n_states, params$n_states))
        for (t in 1:(params$n_obs - 1)) {
          self$omega[t,,] <- matrix(c(0,1,1,0), 2, 2)
        }
      } else {
        self$omega <- array(runif((params$n_obs - 1)*params$n_states*params$n_states), 
                            dim = c(params$n_obs - 1, 
                                    params$n_states, 
                                    params$n_states))
        for (t in 1:(params$n_obs - 1)) {
          diag(self$omega[t,,]) <- 0
          for (k in 1:params$n_states) {
            self$omega[t, k, -k] <- self$omega[t, k, -k]/sum(self$omega[t, k, -k])
          }
        }
      }

      # gamma
      self$gamma <- self$compute_gamma(params, p.array)

    },
    compute_gamma = function(params, p.array) {
      Omega <- self$omega
      K <- params$n_states
      M <- params$max_dwell
      
      # function(Omega, p.array, K, M) {
      # Old code
      n <- nrow(p.array)
      ld <- rep(M, K)
      d <- rep(1:K, ld)
      p <- p.array
      Gamma <- array(0, dim = c(n, sum(ld), sum(ld)))
      for (t in 1:n) {
        for (k in 1:K) {
          pos <- which(d == k)
          if (length(pos) >= 2) {
            if (length(pos) == 2) {
              Gamma[t, min(pos):(max(pos) - 1), (min(pos) + 1):max(pos)] <- 1 - p[t, pos][-length(p[t, pos])]
            } else{
              diag(Gamma[t, min(pos):(max(pos) - 1), (min(pos) + 1):max(pos)]) <- 1 -
                p[t, pos][-length(p[t, pos])]
            }
          }
          Gamma[t, max(pos), max(pos)] <- 1 - p[t, pos][length(p[t, pos])]

          for (i in (1:K)[-k]) {
            Gamma[t, min(pos):max(pos), c(1, which(diff(d) == 1) + 1)[i]] <- p[t, pos] *
              Omega[t, k, i]
          }
        }

      }
      
      self$gamma <- Gamma
    },
    compute_omega = function(params, post.bi.pi.aug) {
      covariates <- params$covariates_omega
      
      alpha_coefs <- list() # initialize a list to save the coefficients of the multinomial model. if K = 2, the list is empty.ci sta 
      if(params$n_states > 2){
        pi_kht <- array(NA, dim = c(params$n_states, params$n_obs-1, params$n_states-1))
        cnt <- 1
        for(k in 1:params$n_states){
          for(t in 1:(params$n_obs-1)){
            cnt2 <- 1
            for(h in 1:params$n_states){
              if(k != h){
                posk <- which(params$state_indices==k)
                posh <- which(params$state_indices==h)
                pi_kht[cnt, t, cnt2] <- sum(post.bi.pi.aug[t, posk, posh])
                cnt2 <- cnt2+1
              }
            }
          }
          
          desing_mat <- as.data.frame(cbind(pi_kht[k,,], covariates[1:(params$n_obs-1)]))
          
          #print(params$n_obs)
          #print(dim(desing_mat))
          #print(params$n_states)
          #print(c(1:params$n_states)[-k])
          
          colnames(desing_mat) <- c(c(1:params$n_states)[-k], "x")
          
          formula_omega <- as.formula("pi ~ .")

          desing_mat <- desing_mat %>%
            pivot_longer(1:(params$n_states-1), names_to = "pi", values_to = "w") %>%
            mutate(pi = factor(pi)) # long format
          
          if(params$n_states == 3){
            # Fit the logistic model
            fit_omega <- suppressWarnings(glm(formula_omega, weights = w, data = desing_mat, family = "binomial"))
            preds <- predict(fit_omega, type = "response") # need the probs!
            preds <- cbind(preds, 1-preds) # add the baseline probability
          } else {
            # Fit the multinomial model
            fit_omega <- suppressWarnings(nnet::multinom(formula_omega, weights = w, data = desing_mat, trace = FALSE))
            preds <- predict(fit_omega, type = "probs") # need the probs!
          }
          
          # Save alpha coefficients of the multinomial model
          alpha_coefs[[k]] <- coef(fit_omega)
          
          # Assign estimated probabilities to Omega
          idextr <- seq(1, nrow(preds), by = params$n_states-1) # index to save: probs repeated due to long format
          appids <- 1:params$n_states
          appids[k] <- 1
          appids[-k] <- (2:params$n_states)
          
          
          self$omega[,k,] <- cbind(0, preds[idextr, ])[, appids]
          cnt <- cnt + 1
        }
      }
    }
  )
)

# ---- End of models_transition_time_varying.R ----

# ---- Start of RcppExports.R ----

# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Gamma_f_cpp <- function(Omega, p_array, K, M, ld, d) {
    .Call(`_hsmmpackage_Gamma_f_cpp`, Omega, p_array, K, M, ld, d)
}

lbackward_forward <- function(lfit, Gamma, Pi, ld, K, M) {
    .Call(`_hsmmpackage_lbackward_forward`, lfit, Gamma, Pi, ld, K, M)
}


# ---- End of RcppExports.R ----

