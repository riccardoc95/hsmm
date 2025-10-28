TorusModel <- R6Class(
  "TorusModel",
  inherit = EmissionModel,

  public = list(
    tau.mu1 = NULL,
    tau.mu2 = NULL,
    tau.k1  = NULL,
    tau.k2  = NULL,
    tau.rho = NULL,
    density = NULL,
    n_states = NULL,
    post.pi  = NULL,

    initialize = function(params, Pi = NULL, post.pi = NULL) {
      super$initialize(params, Pi, post.pi)

      K <- params$n_states
      data <- params$data
      wmat <- post.pi
      if (is.null(wmat)) {
        # fallback uniforme se non hai ancora posteriori
        wmat <- matrix(1/K, nrow = nrow(data), ncol = K)
      }

      self$n_states <- K
      self$post.pi  <- wmat

      # stima iniziale parametri per ogni stato k
      mu1   <- numeric(K)
      mu2   <- numeric(K)
      k1    <- numeric(K)
      k2    <- numeric(K)
      rho12 <- numeric(K)

      for (k in 1:K) {
        w <- wmat[, k]
        sw <- sum(w)
        if (sw < 1e-8) {
          # fallback se nessuna osservazione assegnata
          mu1[k] <- runif(1, 0, 2*pi)
          mu2[k] <- runif(1, 0, 2*pi)
          k1[k]  <- 0.5
          k2[k]  <- 0.5
          rho12[k] <- 0
          next
        }

        # media circolare pesata per ciascuna coordinata angolare
        c1 <- sum(w * cos(data[,1])) / sw
        s1 <- sum(w * sin(data[,1])) / sw
        c2 <- sum(w * cos(data[,2])) / sw
        s2 <- sum(w * sin(data[,2])) / sw

        mu1[k] <- atan2(s1, c1) %% (2*pi)
        mu2[k] <- atan2(s2, c2) %% (2*pi)

        R1 <- sqrt(c1^2 + s1^2)
        R2 <- sqrt(c2^2 + s2^2)

        # "concentrazione proxy", bounded away from 0
        k1[k] <- max(R1, 1e-3)
        k2[k] <- max(R2, 1e-3)

        # correlazione approssimata tra angoli (shiftati sulle loro medie)
        d1 <- data[,1] - mu1[k]
        d2 <- data[,2] - mu2[k]
        # riportiamo in [-pi,pi] perché le differenze angolari sono circolari
        d1 <- (d1 + pi) %% (2*pi) - pi
        d2 <- (d2 + pi) %% (2*pi) - pi

        rho_num <- sum(w * sin(d1) * sin(d2)) / sw
        rho_den <- sqrt( (sum(w * sin(d1)^2)/sw) * (sum(w * sin(d2)^2)/sw) )
        if (!is.finite(rho_den) || rho_den < 1e-8) {
          rho12[k] <- 0
        } else {
          rho12[k] <- max(min(rho_num / rho_den,  0.95), -0.95)
        }
      }

      self$tau.mu1 <- mu1
      self$tau.mu2 <- mu2
      self$tau.k1  <- k1
      self$tau.k2  <- k2
      self$tau.rho <- rho12

      # inizializza densità
      self$compute_density(params)
    },

    # densità torus stabile numericamente
    torus_density = function(theta1, theta2, par) {
      mu1   <- par[1]
      mu2   <- par[2]
      k1    <- max(par[3], 1e-6)
      k2    <- max(par[4], 1e-6)
      rho12 <- max(min(par[5], 0.999), -0.999)

      # costante di normalizzazione "pseudo"
      # ATTENZIONE: nella tua versione usavi Pi come se fosse pi greco.
      const <- (1 - rho12^2) * (1 - k1^2) * (1 - k2^2)
      const <- const / (4 * pi^2)

      # coefficienti
      c0 <- (1 + rho12^2) * (1 + k1^2) * (1 + k2^2) -
        8 * abs(rho12) * k1 * k2

      c1 <- 2 * (1 + rho12^2) * k1 * (1 + k2^2) -
        4 * abs(rho12) * (1 + k1^2) * k2

      c2 <- 2 * (1 + rho12^2) * k2 * (1 + k1^2) -
        4 * abs(rho12) * (1 + k2^2) * k1

      c3 <- -4 * (1 + rho12^2) * k1 * k2 +
        2 * abs(rho12) * (1 + k1^2) * (1 + k2^2)

      c4 <- 2 * rho12 * (1 - k1^2) * (1 - k2^2)

      dtheta1 <- theta1 - mu1
      dtheta2 <- theta2 - mu2

      num <- const
      den <- (
        c0
        - c1 * cos(dtheta1)
        - c2 * cos(dtheta2)
        - c3 * cos(dtheta1) * cos(dtheta2)
        - c4 * (sin(dtheta1) * sin(dtheta2))
      )

      # stabilizzazione numerica:
      den[den < 1e-12] <- 1e-12
      f <- num / den

      # evita 0 o negativo
      f[!is.finite(f) | f <= 0] <- 1e-12
      f
    },

    compute_density = function(params) {
      data <- params$data
      n    <- nrow(data)
      K    <- params$n_states

      dens <- matrix(0, n, K)
      for (k in 1:K) {
        par_k <- c(
          self$tau.mu1[k],
          self$tau.mu2[k],
          self$tau.k1[k],
          self$tau.k2[k],
          self$tau.rho[k]
        )
        dens[, k] <- self$torus_density(data[,1], data[,2], par_k)
      }

      # salva e ritorna
      self$density <- dens
      invisible(self)
    },

    # funzione obiettivo per ottim
    weighted_negloglik = function(par, theta1, theta2, weights) {
      # reparam:
      mu1   <- par[1] %% (2*pi)
      mu2   <- par[2] %% (2*pi)
      k1    <- (tanh(par[3]) + 1)/2
      k2    <- (tanh(par[4]) + 1)/2
      rho12 <- tanh(par[5])

      k1 <- max(k1, 1e-3)
      k2 <- max(k2, 1e-3)
      rho12 <- max(min(rho12, 0.95), -0.95)

      par_phys <- c(mu1, mu2, k1, k2, rho12)

      f <- self$torus_density(theta1, theta2, par_phys)
      val <- -sum(weights * log(f))
      if (!is.finite(val)) val <- .Machine$double.xmax
      val
    },

    update = function(params, Pi = NULL, post.pi = NULL) {
      data <- params$data
      K    <- params$n_states

      if (is.null(post.pi)) {
        post.pi <- self$post.pi
      } else {
        self$post.pi <- post.pi
      }

      for (k in 1:K) {
        w <- post.pi[, k]
        if (sum(w) < 1e-8) next

        init_par <- c(
          self$tau.mu1[k],
          self$tau.mu2[k],
          atanh( 2*self$tau.k1[k] - 1 ),
          atanh( 2*self$tau.k2[k] - 1 ),
          atanh(self$tau.rho[k])
        )

        opt <- try(
          optim(
            par = init_par,
            fn  = self$weighted_negloglik,
            theta1 = data[,1],
            theta2 = data[,2],
            weights = w,
            control = list(maxit = 200)
          ),
          silent = TRUE
        )

        # fallback se optim fallisce
        if (inherits(opt, "try-error") || !is.finite(opt$value)) {
          next
        }

        # ricostruzione parametri "fisici"
        par_hat <- opt$par
        mu1_hat   <- par_hat[1] %% (2*pi)
        mu2_hat   <- par_hat[2] %% (2*pi)
        k1_hat    <- (tanh(par_hat[3]) + 1)/2
        k2_hat    <- (tanh(par_hat[4]) + 1)/2
        rho_hat   <- tanh(par_hat[5])

        k1_hat <- max(k1_hat, 1e-3)
        k2_hat <- max(k2_hat, 1e-3)
        rho_hat <- max(min(rho_hat, 0.95), -0.95)

        self$tau.mu1[k] <- mu1_hat
        self$tau.mu2[k] <- mu2_hat
        self$tau.k1[k]  <- k1_hat
        self$tau.k2[k]  <- k2_hat
        self$tau.rho[k] <- rho_hat
      }

      # aggiorna densità
      self$compute_density(params)
      invisible(self)
    }
  )
)
