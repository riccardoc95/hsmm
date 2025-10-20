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
      for(j in 1:ncol(params$covariates_q)) {
        self$d.covariates[, j] <- rep(
          params$covariates_q[1:(params$n_obs-1), j], 
          each = sum(params$dwell_lengths)
        )
      }
      
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
            drop(t(as.matrix(x_t))%*%as.matrix(beta_x)))
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