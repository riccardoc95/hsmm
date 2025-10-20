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