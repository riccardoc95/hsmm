.as_numeric_matrix <- function(x) {
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

.normalize_probability <- function(x) {
  x[!is.finite(x) | x < 0] <- 0
  if (sum(x) == 0) {
    x[] <- 1 / length(x)
  } else {
    x <- x / sum(x)
  }
  x
}

.normalize_omega <- function(omega) {
  omega <- as.matrix(omega)
  K <- nrow(omega)
  if (ncol(omega) != K) stop("omega must be a square matrix")
  if (K == 1) return(matrix(1, 1, 1))

  omega[!is.finite(omega) | omega < 0] <- 0
  diag(omega) <- 0
  for (k in seq_len(K)) {
    total <- sum(omega[k, ])
    if (total == 0) {
      omega[k, -k] <- 1 / (K - 1)
    } else {
      omega[k, ] <- omega[k, ] / total
    }
  }
  omega
}

.softmax_rows <- function(eta) {
  eta <- as.matrix(eta)
  eta <- eta - apply(eta, 1, max)
  out <- exp(eta)
  out / rowSums(out)
}

.make_initial_post <- function(data, n_states, model_type) {
  n <- nrow(data)
  if (n_states == 1) {
    return(matrix(1, n, 1))
  }

  x <- data
  if (model_type %in% c("poisson", "exponential", "gamma")) {
    x <- log1p(data)
  }
  if (model_type == "torus") {
    x <- do.call(cbind, lapply(seq_len(ncol(data)), function(j) {
      cbind(cos(data[, j]), sin(data[, j]))
    }))
  }

  x <- scale(x)
  x[!is.finite(x)] <- 0

  km <- try(stats::kmeans(x, centers = n_states, nstart = 20), silent = TRUE)
  if (inherits(km, "try-error")) {
    cl <- sample.int(n_states, n, replace = TRUE)
  } else {
    cl <- km$cluster
  }

  post <- matrix(0.05 / (n_states - 1), n, n_states)
  post[cbind(seq_len(n), cl)] <- 0.95
  post
}

.weighted_mean <- function(x, w) {
  sw <- sum(w)
  if (!is.finite(sw) || sw < 1e-10) {
    return(colMeans(x))
  }
  colSums(x * w) / sw
}

.weighted_covariance <- function(x, w, mu = NULL) {
  if (is.null(mu)) {
    mu <- .weighted_mean(x, w)
  }
  sw <- sum(w)
  if (!is.finite(sw) || sw < 1e-10) {
    sw <- nrow(x)
    w <- rep(1, nrow(x))
  }

  centered <- sweep(x, 2, mu, "-")
  sigma <- crossprod(centered * sqrt(w)) / sw
  sigma <- as.matrix(sigma)
  diag(sigma) <- diag(sigma) + 1e-6
  sigma
}

.dmvnorm_basic <- function(x, mean, sigma) {
  x <- as.matrix(x)
  p <- ncol(x)
  sigma <- as.matrix(sigma)
  diag(sigma) <- diag(sigma) + 1e-8

  R <- try(chol(sigma), silent = TRUE)
  if (inherits(R, "try-error")) {
    v <- base::mean(diag(sigma))
    if (!is.finite(v) || v <= 0) v <- 1
    sigma <- diag(p) * v
    R <- chol(sigma)
  }

  z <- sweep(x, 2, mean, "-")
  z <- backsolve(R, t(z), transpose = TRUE)
  quad <- colSums(z^2)
  logdet <- 2 * sum(log(diag(R)))
  exp(-0.5 * (p * log(2 * pi) + logdet + quad))
}

.rmvnorm_basic <- function(n, mean, sigma) {
  mean <- as.numeric(mean)
  p <- length(mean)
  sigma <- as.matrix(sigma)
  diag(sigma) <- diag(sigma) + 1e-8
  R <- try(chol(sigma), silent = TRUE)
  if (inherits(R, "try-error")) {
    R <- chol(diag(p))
  }
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  sweep(z %*% R, 2, mean, "+")
}

.collapse_posterior <- function(post.pi.aug, n_states, max_dwell) {
  n <- nrow(post.pi.aug)
  post.pi <- matrix(0, n, n_states)
  for (k in seq_len(n_states)) {
    pos <- ((k - 1) * max_dwell + 1):(k * max_dwell)
    post.pi[, k] <- rowSums(post.pi.aug[, pos, drop = FALSE])
  }
  post.pi
}

.duration_counts <- function(post.bi.pi.aug, n_states, max_dwell, state) {
  n <- dim(post.bi.pi.aug)[1]
  all_states <- seq_len(n_states * max_dwell)
  from <- ((state - 1) * max_dwell + 1):(state * max_dwell)
  other <- setdiff(all_states, from)

  cases <- numeric(n * max_dwell)
  noncases <- numeric(n * max_dwell)
  id <- 1

  for (t in seq_len(n)) {
    for (m in seq_len(max_dwell)) {
      h <- from[m]
      if (length(other) > 0) {
        cases[id] <- sum(post.bi.pi.aug[t, h, other])
      }
      noncases[id] <- sum(post.bi.pi.aug[t, h, from])
      id <- id + 1
    }
  }

  list(cases = cases, noncases = noncases)
}

.fit_duration_glm <- function(cases, noncases, age, x = NULL, semi = TRUE) {
  keep <- cases + noncases > 1e-12
  if (!any(keep)) {
    return(NULL)
  }

  dat <- data.frame(cases = cases[keep], noncases = noncases[keep])
  terms <- character(0)

  if (semi) {
    dat$age <- age[keep]
    terms <- c(terms, "age")
  }

  if (!is.null(x) && ncol(x) > 0) {
    x <- as.data.frame(x[keep, , drop = FALSE])
    names(x) <- paste0("x", seq_len(ncol(x)))
    dat <- cbind(dat, x)
    terms <- c(terms, names(x))
  }

  rhs <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
  frm <- stats::as.formula(paste("cbind(cases, noncases) ~", rhs))

  fit <- try(
    suppressWarnings(stats::glm(
      frm,
      data = dat,
      family = stats::binomial(link = "cloglog"),
      control = stats::glm.control(maxit = 50)
    )),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  co <- stats::coef(fit)
  co[!is.finite(co)] <- 0
  pmax(pmin(co, 20), -20)
}

.make_lagged_data <- function(y, ar) {
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  rows <- (ar + 1):n
  X <- matrix(1, n - ar, 1 + p * ar)

  for (lag in seq_len(ar)) {
    cols <- (2 + (lag - 1) * p):(1 + lag * p)
    X[, cols] <- y[(ar + 1 - lag):(n - lag), , drop = FALSE]
  }

  list(X = X, Y = y[rows, , drop = FALSE], rows = rows)
}

.weighted_var_fit <- function(X, Y, w, lambda = 0) {
  sw <- sum(w)
  if (!is.finite(sw) || sw < 1e-10) {
    w <- rep(1, nrow(X))
    sw <- nrow(X)
  }

  Xw <- X * sqrt(w)
  Yw <- Y * sqrt(w)
  penalty <- diag(ncol(X)) * lambda
  penalty[1, 1] <- 0

  B <- try(solve(crossprod(Xw) + penalty + diag(ncol(X)) * 1e-8,
                 crossprod(Xw, Yw)), silent = TRUE)
  if (inherits(B, "try-error")) {
    B <- try(qr.solve(Xw, Yw), silent = TRUE)
  }
  if (inherits(B, "try-error")) {
    B <- matrix(0, ncol(X), ncol(Y))
  }

  residuals <- Y - X %*% B
  sigma <- crossprod(residuals * sqrt(w)) / sw
  diag(sigma) <- diag(sigma) + 1e-6

  list(coef = B, sigma = sigma)
}

.predict_softmax_reference <- function(X, B, n_destinations) {
  if (n_destinations == 1) {
    return(matrix(1, nrow(X), 1))
  }
  eta <- cbind(X %*% B, 0)
  .softmax_rows(eta)
}

.fit_transition_softmax <- function(X, Y, start = NULL) {
  J <- ncol(Y)
  P <- ncol(X)
  if (J <= 1 || sum(Y) < 1e-10) {
    return(matrix(0, P, max(J - 1, 0)))
  }

  if (is.null(start) || !all(dim(start) == c(P, J - 1))) {
    start <- matrix(0, P, J - 1)
  }

  objective <- function(par) {
    B <- matrix(par, P, J - 1)
    pr <- .predict_softmax_reference(X, B, J)
    -sum(Y * log(pmax(pr, 1e-12))) + 1e-6 * sum(par^2)
  }

  gradient <- function(par) {
    B <- matrix(par, P, J - 1)
    pr <- .predict_softmax_reference(X, B, J)
    total <- rowSums(Y)
    score <- pr[, seq_len(J - 1), drop = FALSE] * total -
      Y[, seq_len(J - 1), drop = FALSE]
    as.numeric(crossprod(X, score) + 2e-6 * B)
  }

  fit <- try(stats::optim(
    as.numeric(start), objective, gradient,
    method = "BFGS", control = list(maxit = 100, reltol = 1e-8)
  ), silent = TRUE)

  if (inherits(fit, "try-error") || !is.finite(fit$value)) {
    return(start)
  }
  matrix(fit$par, P, J - 1)
}

.transition_counts <- function(post.bi.pi.aug, n_states, max_dwell, from_state) {
  Tm1 <- dim(post.bi.pi.aug)[1]
  from <- ((from_state - 1) * max_dwell + 1):(from_state * max_dwell)
  destinations <- setdiff(seq_len(n_states), from_state)
  out <- matrix(0, Tm1, length(destinations))

  for (j in seq_along(destinations)) {
    state <- destinations[j]
    to <- ((state - 1) * max_dwell + 1):(state * max_dwell)
    for (t in seq_len(Tm1)) {
      out[t, j] <- sum(post.bi.pi.aug[t, from, to])
    }
  }
  out
}

.copy_init_fields <- function(object, values) {
  if (is.null(values) || !is.list(values)) return(invisible(object))
  for (nm in names(values)) {
    ok <- try(object[[nm]] <- values[[nm]], silent = TRUE)
  }
  invisible(object)
}

.rvonmises_basic <- function(n, mu = 0, kappa = 1) {
  if (!is.finite(kappa) || kappa < 1e-8) {
    return(stats::runif(n, 0, 2 * pi))
  }

  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a)) / (2 * kappa)
  r <- (1 + b^2) / (2 * b)
  out <- numeric(n)

  for (i in seq_len(n)) {
    repeat {
      u <- stats::runif(3)
      z <- cos(pi * u[1])
      f <- (1 + r * z) / (r + z)
      c <- kappa * (r - f)
      if (u[2] < c * (2 - c) || log(c / u[2]) + 1 - c >= 0) {
        angle <- if (u[3] > 0.5) acos(f) else -acos(f)
        out[i] <- (mu + angle) %% (2 * pi)
        break
      }
    }
  }
  out
}

.kappa_from_resultant <- function(r) {
  r <- min(max(r, 1e-6), 0.999999)
  if (r < 0.53) {
    k <- 2 * r + r^3 + 5 * r^5 / 6
  } else if (r < 0.85) {
    k <- -0.4 + 1.39 * r + 0.43 / (1 - r)
  } else {
    k <- 1 / (r^3 - 4 * r^2 + 3 * r)
  }
  min(max(k, 1e-3), 100)
}
