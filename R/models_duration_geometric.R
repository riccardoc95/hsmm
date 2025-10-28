#' Geometric Duration Model for HMM/HSMM
#'
#' Implements a **geometric duration model** for Hidden Markov or
#' Hidden Semi-Markov Models (HMM/HSMM).  
#' This model assumes the duration (dwell time) in each hidden state follows
#' a geometric-like distribution derived from a Gompertz hazard base.
#'
#' @description
#' The model estimates per-state intercepts and time effects through
#' logistic regression (GLM) based on expected transitions from the E-step
#' of the EM algorithm.  
#' It is used when no explicit covariates for duration (\code{covariates_q})
#' are supplied.
#'
#' @field p.array Array of per-state and per-duration transition probabilities.
#' @field beta.intercepts Numeric vector of intercept parameters per state.
#' @field beta.time_effects Numeric vector of slope parameters (time effects).
#' @field beta.covariate Duration covariate coefficients (NULL in this model).
#' @field d.state, d.time State and dwell-time indices used in the regression.
#' @field link Link function for the GLM (default: "cloglog").
#'
#' @examples
#' params <- list(n_states = 2, max_dwell = 3, n_obs = 10, semi = TRUE)
#' dm <- DurationModel.Geometric$new(params)
#' str(dm$p.array)
#'
library(R6)

DurationModel.Geometric <- R6Class(
  "DurationModel.Geometric",
  inherit = DurationModel,
  public = list(
    #' @description
    #' Initializes a geometric duration model with random intercepts and time effects.
    #' 
    #' @param params List containing:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states.}
    #'     \item{\code{max_dwell}}{Maximum dwell time per state.}
    #'     \item{\code{n_obs}}{Number of observations.}
    #'   }
    initialize = function(params) {
      super$initialize(params)
      self$beta.covariate <- NULL       # No covariates in this simple model
      self$d.covariates <- NULL
      self$link <- "cloglog"            # Complementary log-log link for hazard
      self$compute_p.array(params)      # Compute initial p.array
    },
    
    #' @description
    #' Computes per-dwell hazard probabilities based on a Gompertz-type rate.
    #' 
    #' @param i Integer vector of dwell-time indices.
    #' @param alpha Numeric intercept parameter.
    #' @param betai Numeric slope (time effect) parameter.
    #' @return Numeric vector of hazard probabilities for each dwell time.
    p.k = function(i, alpha, betai) {
      # Compute Gompertz hazard rate (shifted by +0.5 for stability)
      rate <- self$gomp(i + 0.5, alpha, betai)
      # Convert to probability using exponential form
      p <- 1 - exp(-pmax(rate, .Machine$double.eps))
      p
    },
    
    #' @description
    #' Computes and stores the array of duration probabilities (\code{p.array})
    #' for each state and dwell time across all observations.
    #' 
    #' @param params Model parameter list.
    #' @return Invisibly returns \code{self}.
    compute_p.array = function(params) {
      K <- params$n_states             # Number of hidden states
      M <- params$max_dwell            # Maximum dwell time
      Tm1 <- params$n_obs - 1          # Number of transitions
      
      self$p.array <- array(0, dim = c(Tm1, K * M))
      
      # Compute dwell probabilities for each time and state
      for (t in 1:Tm1) {
        for (k in 1:K) {
          pos <- which(params$state_indices == k)
          self$p.array[t, pos] <- self$p.k(
            1:length(pos),
            alpha = self$beta.intercepts[k],
            betai = self$beta.time_effects[k]
          )
        }
      }
    },
    
    #' @description
    #' Estimates duration model parameters (\code{beta.intercepts}, \code{beta.time_effects})
    #' via logistic regression using posterior transition counts from the E-step.
    #' 
    #' @param params Model parameter list.
    #' @param post.bi.pi.aug 3D array of posterior transition probabilities
    #'   between augmented states (time × from × to).
    #' @return Invisibly returns \code{self}.
    compute_beta = function(params, post.bi.pi.aug) {
      n <- params$n_obs - 1
      d <- rep(1:params$n_states, each = params$max_dwell)
      
      cases_fin <- noncases_fin <- numeric()
      
      # ---- Count transitions and non-transitions ----
      for (t in 1:n) {
        for (h in 1:params$n_states) {
          posh <- which(d == h)
          others <- which(d != h)
          
          postbipi_1 <- post.bi.pi.aug[t, posh, others, drop = FALSE]  # exits
          postbipi_2 <- post.bi.pi.aug[t, posh, posh, drop = FALSE]    # stays
          
          cases_fin <- c(cases_fin, rowSums(postbipi_1))
          noncases_fin <- c(noncases_fin, rowSums(postbipi_2))
        }
      }
      
      outno0 <- (cases_fin + noncases_fin) != 0
      if (!any(outno0)) return(invisible(self))
      
      # ---- Build regression dataset ----
      d_matrix <- data.frame(
        state = rep(self$d.state, times = params$n_obs - 1),
        time = rep(self$d.time, times = params$n_obs - 1)
      )
      if (!is.null(self$d.covariates)) {
        d_matrix <- cbind(
          d_matrix,
          setNames(as.data.frame(self$d.covariates), paste0("V", 1:ncol(self$d.covariates)))
        )
      }
      
      df_for_glm <- data.frame(
        cases = cases_fin[outno0],
        noncases = noncases_fin[outno0],
        d_matrix[outno0, , drop = FALSE]
      )
      
      # ---- Fit GLM model ----
      formula_dwell <- as.formula("cbind(cases, noncases) ~ factor(state) * time")
      fit_beta <- suppressWarnings(
        glm(formula_dwell,
            family = binomial(link = self$link),
            data = df_for_glm)
      )
      
      coefs <- coef(fit_beta)
      if (is.null(coefs)) return(invisible(self))
      
      # ---- Extract intercepts per state ----
      beta_intercepts <- numeric(params$n_states)
      beta_intercepts[1] <- coefs["(Intercept)"]
      for (i in 2:params$n_states) {
        nm <- paste0("factor(state)", i)
        beta_intercepts[i] <- beta_intercepts[1] + ifelse(!is.na(coefs[nm]), coefs[nm], 0)
      }
      
      # ---- Extract time effects ----
      time_id <- grep("time", names(coefs))
      beta_time_effects <- numeric(params$n_states)
      beta_time_effects[] <- if (length(time_id) > 0) coefs[time_id][1] else 0
      
      # ---- Store updated parameters ----
      self$beta.intercepts <- beta_intercepts
      self$beta.time_effects <- if (params$semi) beta_time_effects else rep(0, length(beta_intercepts))
      self$beta.covariate <- NULL
    }
  )
)
