\name{simulate_hsmm}
\alias{simulate_hsmm}
\alias{simulate.hsmm_fit}
\title{Simulate a hidden semi-Markov series}
\usage{
simulate_hsmm(
  n = 300,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 1,
  semi = TRUE,
  max_dwell = 20,
  covariates_omega = NULL,
  covariates_q = NULL,
  Pi = NULL,
  omega = NULL,
  duration_intercept = NULL,
  duration_age_effect = NULL,
  duration_coef = NULL,
  transition_coef = NULL,
  emission = NULL,
  ar = 1,
  seed = NULL,
  family = NULL
)

\method{simulate}{hsmm_fit}(
  object,
  nsim = 1,
  seed = NULL,
  n = object$input_params$n_obs,
  covariates_omega = NULL,
  covariates_q = NULL,
  ...
)
}
\arguments{
\item{n}{Number of observations.}
\item{n_states}{Number of latent states.}
\item{model_type}{Emission model name.}
\item{family}{Old alias for \code{model_type}.}
\item{n_dim}{Number of observed variables.}
\item{semi}{Use an age-dependent duration hazard.}
\item{max_dwell}{Number of duration states; the last represents the tail.}
\item{covariates_omega}{Optional transition covariates.}
\item{covariates_q}{Optional duration covariates.}
\item{Pi}{Initial state probabilities.}
\item{omega}{Conditional transition matrix.}
\item{duration_intercept}{Cloglog duration intercepts.}
\item{duration_age_effect}{Duration age effects.}
\item{duration_coef}{Duration covariate coefficients.}
\item{transition_coef}{List of transition coefficient matrices.}
\item{emission}{List of emission parameters.}
\item{ar}{VAR order.}
\item{seed}{Optional random seed.}
\item{object}{A fitted \code{hsmm_fit} object.}
\item{nsim}{Number of series to simulate.}
\item{...}{Unused.}
}
\value{A list of class \code{hsmm_simulation}, or a list of simulations.}
\description{
Simulates latent states and observations from supplied or fitted parameters.
}
\examples{
sim <- simulate_hsmm(n = 100, n_states = 3, seed = 2)
head(sim$data)
head(sim$state)
}
