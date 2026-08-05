\name{fit_hsmm}
\alias{fit_hsmm}
\title{Fit a hidden semi-Markov model}
\usage{
fit_hsmm(
  data,
  n_states,
  covariates_omega = NULL,
  covariates_q = NULL,
  semi = TRUE,
  max_dwell = NULL,
  model_type = "gaussian",
  max_iter = 100,
  tol = 1e-05,
  verbose = TRUE,
  init = NULL,
  seed = NULL,
  ar = 1,
  lambda = 0,
  family = NULL
)
}
\arguments{
\item{data}{Numeric matrix or data frame. Rows are time points.}
\item{n_states}{Number of latent states.}
\item{covariates_omega}{Optional transition covariates.}
\item{covariates_q}{Optional duration covariates.}
\item{semi}{Use an HSMM. If false, \code{max_dwell} is one.}
\item{max_dwell}{Number of augmented duration states per latent state.}
\item{model_type}{Emission model name.}
\item{family}{Old alias for \code{model_type}.}
\item{max_iter}{Maximum EM iterations.}
\item{tol}{Convergence tolerance.}
\item{verbose}{Print iteration information.}
\item{init}{Optional list of initial values.}
\item{seed}{Optional random seed.}
\item{ar}{VAR order for \code{vargaussian}.}
\item{lambda}{Ridge penalty for \code{vargaussian}.}
}
\value{An object of class \code{hsmm_fit}.}
\description{
Fits a hidden semi-Markov model by EM using an augmented state space.
}
\examples{
sim <- simulate_hsmm(n = 80, n_states = 2, seed = 1)
fit <- fit_hsmm(sim$data, 2, max_dwell = 8, max_iter = 5,
                verbose = FALSE, seed = 1)
fit
}
