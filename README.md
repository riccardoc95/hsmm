# hsmm

A small R package for fitting and simulating hidden semi-Markov models (HSMMs)
and hidden Markov models (HMMs). The code is intentionally direct: R6 objects
hold the duration, transition and emission pieces, while the augmented-state
forward-backward algorithm is written in C++ with RcppArmadillo.



## Installation

From the repository directory:

```r
install.packages(c("R6", "Rcpp", "RcppArmadillo"))
install.packages(".", repos = NULL, type = "source")
```

For development:

```r
install.packages("testthat")
testthat::test_local()
source(system.file("tests", "manual", "run_all_tests.R", package = "hsmm"))
```

## First check: simulate and fit

```r
library(hsmm)

sim <- simulate_hsmm(
  n = 250,
  n_states = 3,
  model_type = "gaussian",
  n_dim = 2,
  max_dwell = 12,
  seed = 42
)

fit <- fit_hsmm(
  data = sim$data,
  n_states = 3,
  model_type = "gaussian",
  max_dwell = 12,
  max_iter = 50,
  seed = 42
)

fit
summary(fit)
plot(fit)
head(fit$posterior)
head(fit$decoded_state)
```

## Covariates

```r
set.seed(1)
n <- 300
x_q <- matrix(rnorm(n), ncol = 1)
x_omega <- cbind(sin(seq(0, 6 * pi, length.out = n)))
transition_coef <- list(
  matrix(c(0, 1.0), 2, 1),
  matrix(c(0, -0.8), 2, 1),
  matrix(c(0, 0.5), 2, 1)
)

sim <- simulate_hsmm(
  n = n,
  n_states = 3,
  covariates_q = x_q,
  covariates_omega = x_omega,
  duration_coef = matrix(c(0.8, -0.4, 0.2), 3, 1),
  transition_coef = transition_coef,
  seed = 1
)

fit <- fit_hsmm(
  sim$data,
  n_states = 3,
  covariates_q = x_q,
  covariates_omega = x_omega,
  max_dwell = 15,
  max_iter = 30,
  verbose = FALSE,
  seed = 1
)
```

## Emission models

The argument `model_type` accepts:

- `gaussian`
- `poisson`
- `exponential`
- `gamma`
- `beta`
- `student`
- `torus`
- `vargaussian` (TODO!)

The old argument name `family` is kept as an alias for `model_type`.

## Main objects returned by `fit_hsmm`

- `posterior`: posterior probability of each latent state;
- `decoded_state`: maximum-posterior state sequence;
- `loglik` and `final_loglik`;
- `duration.model`, `transition.model`, `emission.model`;
- `posterior_augmented` and `pairwise_posterior` for lower-level work.

## References

Lagona, F. and Mingione, M. (2025). *Nonhomogeneous Hidden Semi-Markov
Models for Toroidal Data*. Journal of the Royal Statistical Society Series C,
74(1), 142-166. https://doi.org/10.1093/jrsssc/qlae049

Mingione, M., Di Loro, P. A., Lagona, F. and Maruotti, A. (2025).
*Environmental Risk Assessment via Nonhomogeneous Hidden Semi-Markov Models
with Penalized Vector Auto-Regression*. arXiv:2509.14387.
