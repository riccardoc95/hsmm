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
      
      #Gamma_f_cpp(Omega, p.array, K, M, ld, d)
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