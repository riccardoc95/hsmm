#' Backward-Forward Algorithm for Hidden Markov Models (HMM)
#'
#' This function performs the forward-backward (also known as backward-forward) algorithm
#' for a Hidden Markov Model (HMM). It computes the forward (`alpha`) and backward (`Beta`)
#' probabilities, as well as the smoothed posterior probabilities of the hidden states.
#'
#' @param fit A numeric matrix (n x S) containing the emission probabilities for each time step and state.
#' @param Gamma A 3-dimensional array (n-1 x S x S) of state transition probabilities between time steps.
#' @param Pi A numeric vector of length S with the initial state probabilities.
#' @param K An integer, the number of latent groups (or sub-models).
#' @param M An integer, the number of states per group.
#'
#' @return A list containing:
#' \describe{
#'   \item{l}{The log-likelihood of the observed data under the model.}
#'   \item{post.pi.aug}{A matrix (n x S) of posterior probabilities for each state at each time step.}
#'   \item{post.bi.pi.aug}{A 3D array (n-1 x S x S) of pairwise posterior probabilities for transitions.}
#' }
#'
#' @details
#' This implementation includes scaling at each time step (using the vector `C`) to avoid numerical underflow.
#' The function iteratively computes forward and backward probabilities and combines them to obtain smoothed
#' posterior probabilities for both states and transitions.
#'
#' @examples
#' # Example usage:
#' # backward.forward(fit = matrix(runif(100), 10, 10),
#' #                  Gamma = array(runif(900), c(9, 10, 10)),
#' #                  Pi = rep(1/10, 10),
#' #                  K = 2, M = 5)
#'
backward.forward <- function(fit, Gamma, Pi, K, M) {

  out <- backward_forward(fit, Gamma, Pi, K, M)
  return(list(l = out$l,
              post.pi.aug = out$post.pi.aug,
              post.bi.pi.aug = out$post.bi.pi.aug))


  # cat("dim(fit) =", paste(dim(fit), collapse="x"), "\n")
  # cat("dim(Gamma) =", paste(dim(Gamma), collapse="x"), "\n")
  # cat("length(Pi) =", length(Pi), "\n")
  # cat("K =", K, "\n")
  # cat("M =", M, "\n")
  # cat("K*M =", K*M, "\n")
  #
  #
  # # Compute number of states per group
  # ld <- rep(M, K)
  # n <- nrow(fit)  # Number of observations (time steps)
  # S <- sum(ld)    # Total number of states
  #
  # # Initialize matrices
  # psi <- alpha <- Beta <- matrix(0, n, S)  # Forward and backward probabilities
  # C <- numeric(n)                          # Scaling coefficients
  #
  # # ---- Forward pass (alpha) ----
  # psi[1, ] <- Pi                           # Initial state probabilities
  # psi_per_fit <- psi[1, ] * fit[1, ]       # Weighted by emission probabilities
  # C[1] <- sum(psi_per_fit)                 # Scaling factor for step 1
  #
  # # Prevent numerical underflow
  # if (C[1] < .Machine$double.eps) C[1] <- .Machine$double.eps
  #
  # # Normalize the first forward probabilities
  # alpha[1, ] <- psi_per_fit / C[1]
  #
  # # Iterate over all time steps
  # for (i in 2:n) {
  #   # Predict next state probabilities using transition matrix
  #   psi[i, ] <- as.numeric(alpha[i-1, ] %*% Gamma[i-1, , ])
  #
  #   # Multiply by emission probabilities for current step
  #   psi_per_fit <- psi[i, ] * fit[i, ]
  #
  #   # Compute scaling factor and prevent underflow
  #   C[i] <- sum(psi_per_fit)
  #   if (C[i] < .Machine$double.eps) C[i] <- .Machine$double.eps
  #
  #   # Normalize forward probabilities
  #   alpha[i, ] <- psi_per_fit / C[i]
  # }
  #
  # # Compute log-likelihood as sum of log scaling factors
  # l <- sum(log(C))
  #
  # # ---- Backward pass (Beta) ----
  # Beta[n, ] <- 1  # Initialize last backward probability to 1
  #
  # # Iterate backward from n-1 to 1
  # for (i in (n-1):1) {
  #   # Compute backward probability for each state
  #   Beta[i, ] <- (Gamma[i, , ] %*% (fit[i+1, ] * Beta[i+1, ])) / C[i+1]
  # }
  #
  # # ---- Posterior probabilities ----
  # alphaBeta <- alpha * Beta
  # post.pi.aug <- alphaBeta / rowSums(alphaBeta)  # State posterior
  #
  # # ---- Pairwise posterior probabilities ----
  # post.bi.pi.aug <- array(0, dim = c(n-1, S, S))
  #
  # # Compute posterior for all state transitions
  # for (h in 1:S) {
  #   for (k in 1:S) {
  #     post.bi.pi.aug[, h, k] <- alpha[1:(n-1), h] * Gamma[, h, k] *
  #       fit[2:n, k] * Beta[2:n, k] / C[2:n]
  #   }
  # }
  #
  # # Return log-likelihood and posterior probabilities
  # return(list(l = l, post.pi.aug = post.pi.aug, post.bi.pi.aug = post.bi.pi.aug))
}
