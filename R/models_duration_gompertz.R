#' Gompertz Duration Model for HSMM
#'
#' Implements a **covariate-dependent Gompertz duration model** for Hidden
#' Semi-Markov Models (HSMMs).  
#' This class extends the base `DurationModel` and introduces covariates
#' that influence the duration (dwell-time) distribution in each state
#' through a complementary log-log link (cloglog).
#'
#' @description
#' The Gompertz model is used when duration covariates
#' (\code{covariates_q}) are provided in the parameter list.
#' It models the hazard rate (probability of exiting a state)
#' as a function of time and covariates using a Gompertz-type hazard base.
#'
#' The parameters are estimated iteratively within the EM algorithm:
#' - **E-step**: posterior transitions are computed.
#' - **M-step**: parameters (\code{beta.intercepts}, \code{beta.time_effects},
#'   and \code{beta.covariate}) are updated via GLM fitting.
#'
#' @field p.array 3D array of per-state, per-duration, per-time exit probabilities.
#' @field beta.intercepts Vector of per-state intercepts.
#' @field beta.time_effects Vector of per-state time effects.
#' @field beta.covariate Matrix of per-state covariate coefficients.
#' @field d.state, d.time Internal dwell-time indexing vectors.
#' @field link Link function used in GLM (always "cloglog").
#'
#' @examples
#' params <- list(
#'   n_states = 2,
#'   max_dwell = 3,
#'   n_obs = 10,
#'   semi = TRUE,
#'   covariates_q = matrix(rnorm(10 * 2), ncol = 2)
#' )
#' dm <- DurationModel.Gompertz$new(params)
#' str(dm$p.array)
#'
library(R6)

DurationModel.Gompertz <- R6Class(
  "DurationModel.Gompertz",
  inherit = DurationModel,
  public = list(
    #' @description
    #' Initializes the Gompertz duration model with random coefficients and covariates.
    #' 
    #' @param params List containing:
    #'   \itemize{
    #'     \item{\code{n_states}}{Number of hidden states.}
    #'     \item{\code{max_dwell}}{Maximum dwell time.}
    #'     \item{\code{n_obs}}{Number of observations.}
    #'     \item{\code{covariates_q}}{Matrix of dwell-time covariates.}
    #'   }
    initialize = function(params) {
      super$initialize(params)
      
      K <- params$n_states             # Number of hidden states
      Q <- ncol(params$covariates_q)   # Number of covariates
      
      # Initialize per-state covariate coefficients
      self$beta.covariate <- matrix(runif(K * Q, -0.5, 0.5), nrow = K, ncol = Q)
      
      # Repeat covariates for each dwell-time instance across all states
      self$d.covariates <- matrix(
        rep(as.matrix(params$covariates_q[1:(params$n_obs - 1), ]),
            each = K * params$max_dwell / K),
        ncol = Q
      )
      
      self$link <- "cloglog"           # Complementary log-log link
      self$compute_p.array(params)     # Compute initial probability array
    },
    
    #' @description
    #' Computes the per-time hazard probability for the Gompertz model.
    #' 
    #' @param i Dwell-time indices.
    #' @param alpha State-specific intercept.
    #' @param betai Time effect coefficient.
    #' @param beta_x Vector of covariate coefficients for the state.
    #' @param x_t Covariate vector at time t.
    #' @return Numeric vector of hazard probabilities.
    p.k = function(i, alpha, betai, beta_x, x_t) {
      # Linear predictor from covariates
      eta <- drop(t(as.matrix(x_t)) %*% as.matrix(beta_x))
      
      # Apply exponentiation, clamped for numerical stability
      exp_term <- exp(pmin(eta, 700))
      
      # Compute Gompertz hazard rate
      rate <- self$gomp(i + 0.5, alpha, betai)
      
      # Convert to exit probability
      1 - exp(-pmax(rate, .Machine$double.eps) * exp_term)
    },
    
    #' @description
    #' Computes the array of dwell-time probabilities across states and times,
    #' combining Gompertz hazard dynamics and covariate effects.
    #' 
    #' @param params Model parameters.
    #' @return Invisibly returns \code{self}.
    compute_p.array = function(params) {
      K <- params$n_states
      M <- params$max_dwell
      Tm1 <- params$n_obs - 1
      
      self$p.array <- array(0, dim = c(Tm1, K * M))
      
      # Compute p.array entries for each state and time
      for (t in 1:Tm1) {
        for (k in 1:K) {
          pos <- which(params$state_indices == k)
          self$p.array[t, pos] <- self$p.k(
            1:length(pos),
            alpha = self$beta.intercepts[k],
            betai = self$beta.time_effects[k],
            beta_x = self$beta.covariate[k, ],
            x_t = params$covariates_q[t, ]
          )
        }
      }
    },
    
    #' @description
    #' Updates the parameters of the Gompertz duration model via GLM fitting.
    #' 
    #' @param params Model parameter list.
    #' @param post.bi.pi.aug Posterior bi-state transition probabilities
    #'   (from the E-step of the EM algorithm).
    #' @return Invisibly returns \code{self}.
    compute_beta = function(params, post.bi.pi.aug) {
      n <- params$n_obs - 1
      K <- params$n_states
      M <- params$max_dwell
      Q <- ncol(params$covariates_q)
      
      d <- rep(1:K, each = M)
      cases_fin <- noncases_fin <- numeric()
      
      # ---- Collect exit and stay counts ----
      for (t in 1:n) {
        for (h in 1:K) {
          posh <- which(d == h)
          others <- which(d != h)
          postbipi_1 <- post.bi.pi.aug[t, posh, others, drop = FALSE]  # exits
          postbipi_2 <- post.bi.pi.aug[t, posh, posh, drop = FALSE]    # stays
          cases_fin <- c(cases_fin, rowSums(postbipi_1))
          noncases_fin <- c(noncases_fin, rowSums(postbipi_2))
        }
      }
      
      # Skip if all counts are zero
      outno0 <- (cases_fin + noncases_fin) != 0
      if (!any(outno0)) return(invisible(self))
      
      cov_names <- paste0("V", 1:Q)
      
      # ---- Build data frame for regression ----
      d_matrix <- data.frame(
        state = rep(self$d.state, times = params$n_obs - 1),
        time = rep(self$d.time, times = params$n_obs - 1)
      )
      if (!is.null(self$d.covariates)) {
        d_matrix <- cbind(
          d_matrix,
          setNames(as.data.frame(self$d.covariates),
                   paste0("V", 1:ncol(self$d.covariates)))
        )
      }
      
      df_for_glm <- data.frame(
        cases = cases_fin[outno0],
        noncases = noncases_fin[outno0],
        d_matrix[outno0, , drop = FALSE]
      )
      
      # ---- Construct GLM formula dynamically ----
      formula_str <- paste0(
        "cbind(cases, noncases) ~ factor(state)*time + ",
        paste0("factor(state)*", cov_names, collapse = " + ")
      )
      formula_dwell <- as.formula(formula_str)
      
      # ---- Fit GLM model ----
      fit_beta <- suppressWarnings(
        glm(formula_dwell,
            family = binomial(link = self$link),
            data = df_for_glm)
      )
      
      coefs <- coef(fit_beta)
      if (is.null(coefs)) return(invisible(self))
      
      # ---- Extract intercepts ----
      beta_intercepts <- numeric(K)
      beta_intercepts[1] <- coefs["(Intercept)"]
      for (i in 2:K) {
        nm <- paste0("factor(state)", i)
        beta_intercepts[i] <- beta_intercepts[1] + ifelse(!is.na(coefs[nm]), coefs[nm], 0)
      }
      
      # ---- Extract time effects ----
      time_id <- grep("time", names(coefs))
      beta_time_effects <- numeric(K)
      beta_time_effects[] <- if (length(time_id) > 0) coefs[time_id][1] else 0
      
      # ---- Extract covariate effects ----
      beta_covar <- matrix(0, nrow = K, ncol = Q)
      for (j in 1:Q) {
        base_name <- paste0("V", j)
        base_coef <- coefs[grep(paste0("\\b", base_name, "\\b"), names(coefs))][1]
        if (is.na(base_coef)) base_coef <- 0
        beta_covar[1, j] <- base_coef
        for (k in 2:K) {
          nm <- paste0("factor(state)", k, ":", base_name)
          add_coef <- coefs[nm]
          if (is.na(add_coef)) add_coef <- 0
          beta_covar[k, j] <- beta_covar[1, j] + add_coef
        }
      }
      
      # ---- Update model parameters ----
      self$beta.intercepts <- beta_intercepts
      self$beta.time_effects <- if (params$semi) beta_time_effects else rep(0, length(beta_intercepts))
      self$beta.covariate <- beta_covar
    }
  )
)
