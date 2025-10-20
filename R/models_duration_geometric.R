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