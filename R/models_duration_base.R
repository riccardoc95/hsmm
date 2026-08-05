# Base class for the dwell-time part of the model.
DurationModel <- R6::R6Class(
  "DurationModel",
  public = list(
    p.array = NULL,
    beta.intercepts = NULL,
    beta.time_effects = NULL,
    beta.covariate = NULL,
    link = "cloglog",

    initialize = function(params) {
      K <- params$n_states
      Q <- if (is.null(params$covariates_q)) 0 else ncol(params$covariates_q)

      self$beta.intercepts <- stats::runif(K, -2.5, -1)
      self$beta.time_effects <- if (params$semi) {
        stats::runif(K, 0.02, 0.15)
      } else {
        rep(0, K)
      }
      self$beta.covariate <- matrix(0, K, Q)

      if (!is.null(params$init$duration)) {
        .copy_init_fields(self, params$init$duration)
      }
      self$compute_p.array(params)
    },

    compute_p.array = function(params) {
      K <- params$n_states
      M <- params$max_dwell
      Tm1 <- params$n_obs - 1
      Q <- ncol(self$beta.covariate)
      x <- params$covariates_q

      self$p.array <- matrix(0, Tm1, K * M)
      for (k in seq_len(K)) {
        for (m in seq_len(M)) {
          eta <- rep(self$beta.intercepts[k] +
                       self$beta.time_effects[k] * (m + 0.5), Tm1)
          if (Q > 0) {
            eta <- eta + drop(x[seq_len(Tm1), , drop = FALSE] %*%
                                self$beta.covariate[k, ])
          }
          p <- 1 - exp(-exp(pmax(pmin(eta, 20), -20)))
          self$p.array[, (k - 1) * M + m] <- pmin(pmax(p, 1e-6), 1 - 1e-6)
        }
      }
      invisible(self)
    },

    compute_beta = function(params, post.bi.pi.aug) {
      K <- params$n_states
      M <- params$max_dwell
      Tm1 <- params$n_obs - 1
      Q <- ncol(self$beta.covariate)

      if (K == 1) return(invisible(self))

      age <- rep(seq_len(M) + 0.5, times = Tm1)
      x <- NULL
      if (Q > 0) {
        ids <- rep(seq_len(Tm1), each = M)
        x <- params$covariates_q[ids, , drop = FALSE]
      }

      for (k in seq_len(K)) {
        counts <- .duration_counts(post.bi.pi.aug, K, M, k)
        co <- .fit_duration_glm(
          counts$cases, counts$noncases, age, x, params$semi
        )
        if (is.null(co)) next

        self$beta.intercepts[k] <- unname(co["(Intercept)"])
        if (!is.finite(self$beta.intercepts[k])) self$beta.intercepts[k] <- -2

        if (params$semi && "age" %in% names(co)) {
          self$beta.time_effects[k] <- unname(co["age"])
        } else {
          self$beta.time_effects[k] <- 0
        }

        if (Q > 0) {
          for (j in seq_len(Q)) {
            nm <- paste0("x", j)
            self$beta.covariate[k, j] <- if (nm %in% names(co)) co[nm] else 0
          }
        }
      }
      invisible(self)
    }
  )
)
