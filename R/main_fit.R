#' Fit a hidden semi-Markov model
#'
#' Fits a small HSMM/HMM by EM. Duration hazards use a complementary log-log
#' regression and transition covariates use conditional multinomial logits.
#'
#' @param data Numeric matrix or data frame. Rows are time points.
#' @param n_states Number of latent states.
#' @param covariates_omega Optional transition covariates.
#' @param covariates_q Optional duration covariates.
#' @param semi If `FALSE`, `max_dwell` is forced to one.
#' @param max_dwell Number of augmented duration states per latent state.
#' @param model_type Emission model: `gaussian`, `poisson`, `exponential`,
#'   `gamma`, `beta`, `student`, `torus`, or `vargaussian`.
#' @param family Old alias for `model_type`.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance on the log-likelihood difference.
#' @param verbose Print one line per iteration.
#' @param init Optional list of starting values.
#' @param seed Optional random seed.
#' @param ar Autoregressive order for `vargaussian`.
#' @param lambda Ridge penalty for `vargaussian`.
#'
#' @return An object of class `hsmm_fit`.
#' @export
fit_hsmm <- function(
    data,
    n_states,
    covariates_omega = NULL,
    covariates_q = NULL,
    semi = TRUE,
    max_dwell = NULL,
    model_type = "gaussian",
    max_iter = 100,
    tol = 1e-5,
    verbose = TRUE,
    init = NULL,
    seed = NULL,
    ar = 1,
    lambda = 0,
    family = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (!is.null(family)) model_type <- family
  model_type <- tolower(model_type)
  if (model_type %in% c("studentt", "tstudent")) model_type <- "student"
  if (model_type %in% c("var", "var_gaussian")) model_type <- "vargaussian"
  allowed_models <- c("gaussian", "poisson", "exponential", "gamma",
                      "beta", "student", "torus", "vargaussian")
  if (!model_type %in% allowed_models) stop("unknown model_type")

  data <- .as_numeric_matrix(data)
  n_obs <- nrow(data)
  n_states <- as.integer(n_states)
  ar <- as.integer(ar)

  if (n_obs < 2) stop("data must contain at least two rows")
  if (n_states < 1 || n_states > n_obs) stop("n_states is not valid")
  if (any(!is.finite(data))) stop("data contains non-finite values")
  if (max_iter < 1) stop("max_iter must be positive")
  if (!is.finite(tol) || tol < 0) stop("tol is not valid")
  if (!is.finite(lambda) || lambda < 0) stop("lambda is not valid")

  if (!semi) {
    max_dwell <- 1L
  } else if (is.null(max_dwell)) {
    max_dwell <- min(20L, max(2L, as.integer(n_obs / 10)))
  }
  max_dwell <- as.integer(max_dwell)
  if (max_dwell < 1) stop("max_dwell must be positive")

  prepare_covariates <- function(x, name) {
    if (is.null(x)) return(matrix(numeric(0), n_obs, 0))
    x <- .as_numeric_matrix(x)
    if (nrow(x) == n_obs - 1) x <- rbind(x, x[nrow(x), , drop = FALSE])
    if (nrow(x) != n_obs) stop(name, " must have nrow(data) rows")
    if (any(!is.finite(x))) stop(name, " contains non-finite values")
    x
  }

  covariates_omega <- prepare_covariates(covariates_omega, "covariates_omega")
  covariates_q <- prepare_covariates(covariates_q, "covariates_q")

  if (model_type == "poisson" &&
      any(data < 0 | abs(data - round(data)) > 1e-8)) {
    stop("poisson data must be non-negative integers")
  }
  if (model_type %in% c("exponential", "gamma") && any(data <= 0)) {
    stop(model_type, " data must be positive")
  }
  if (model_type == "beta" && any(data <= 0 | data >= 1)) {
    stop("beta data must be inside (0, 1)")
  }
  if (model_type == "vargaussian" && (ar < 1 || n_obs <= ar + 1)) {
    stop("ar is too large for the supplied data")
  }

  if (is.null(init)) init <- list()
  initial_post <- init$post.pi
  if (is.null(initial_post)) {
    initial_post <- .make_initial_post(data, n_states, model_type)
  } else {
    initial_post <- as.matrix(initial_post)
    if (!all(dim(initial_post) == c(n_obs, n_states))) {
      stop("init$post.pi must have nrow(data) rows and n_states columns")
    }
    initial_post[!is.finite(initial_post) | initial_post < 0] <- 0
    initial_post <- initial_post / rowSums(initial_post)
    initial_post[!is.finite(initial_post)] <- 1 / n_states
  }

  params <- list(
    data = data,
    n_states = n_states,
    n_obs = n_obs,
    covariates_omega = covariates_omega,
    covariates_q = covariates_q,
    semi = isTRUE(semi),
    max_dwell = max_dwell,
    dwell_lengths = rep(max_dwell, n_states),
    state_indices = rep(seq_len(n_states), each = max_dwell),
    model_type = model_type,
    max_iter = as.integer(max_iter),
    tol = tol,
    verbose = isTRUE(verbose),
    init = init,
    initial_post = initial_post,
    seed = seed,
    ar = ar,
    lambda = lambda
  )

  duration.model <- DurationModelFactory$new()$create(params)
  transition.model <- TransitionModelFactory$new()$create(
    params, duration.model$p.array
  )
  emission.model <- EmissionModelFactory$new()$create(
    params, transition.model$Pi, transition.model$post.pi
  )

  em <- em.algorithm(params, emission.model, transition.model, duration.model)
  decoded <- max.col(em$posterior, ties.method = "first")

  out <- list(
    call = match.call(),
    input_params = params,
    emission.model = emission.model,
    transition.model = transition.model,
    duration.model = duration.model,
    loglik = em$loglik,
    final_loglik = em$final_loglik,
    posterior = em$posterior,
    posterior_augmented = em$posterior_augmented,
    pairwise_posterior = em$pairwise_posterior,
    decoded_state = decoded,
    converged = em$converged,
    iterations = em$iterations
  )
  class(out) <- "hsmm_fit"
  out
}
