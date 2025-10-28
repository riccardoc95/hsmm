# hsmm: Hidden Semi-Markov Models with Covariates in R

> Framework for fitting non-homogeneous Hidden Semi-Markov Models (HSMMs) and Hidden Markov Models (HMMs) with covariates, explicit duration modeling, and flexible emission distributions.

---

## Overview

`hsmm` is an R package for estimating **Hidden Semi-Markov Models (HSMMs)** and **Hidden Markov Models (HMMs)** under a unified and extensible framework.  
It supports **time-varying covariates**, **explicit dwell-time distributions**, and **various emission models** (Gaussian, Poisson, Gamma, Torus, etc.), allowing for both homogeneous and non-homogeneous latent dynamics.

This package implements the methodology introduced in:

- Lagona & Mingione (2025). *Nonhomogeneous Hidden Semi-Markov Models for Toroidal Data*.  
  *Journal of the Royal Statistical Society: Series C (Applied Statistics)*, 74(1), 142–166.  
  [https://doi.org/10.1093/jrsssc/qlae049](https://doi.org/10.1093/jrsssc/qlae049)

- Mingione, Di Loro, Lagona & Maruotti (2025).  
  *Environmental Risk Assessment via Nonhomogeneous Hidden Semi-Markov Models with Penalized Vector Auto-Regression*.  
  *arXiv preprint*, [arXiv:2509.14387v1](https://arxiv.org/abs/2509.14387v1).

---


## Installation

```r
# Install dependencies
install.packages(c("Rcpp", "R6", "stats", "MASS", "mvtnorm", "Matrix"))

# Install from GitHub (after building)
devtools::install_github("yourusername/hsmm")
````

---

## Basic Usage

```r
library(hsmm)

# Example: Fit a Gaussian HSMM with 3 hidden states
model <- fit_hsmm(
  data = your_timeseries,
  n_states = 3,
  semi = TRUE,
  family = "gaussian",
  covariates_omega = weather_covariates,
  covariates_q = wind_speed
)

summary(model)
plot(model)
```

---

## Architecture

The package is structured into modular components:

| Component             | Description                                | Key File                                            |
| --------------------- | ------------------------------------------ | --------------------------------------------------- |
| **Emission Models**   | Define how observations are generated      | `R/models_emission_*.R`                             |
| **Duration Models**   | Define dwell-time distributions            | `R/models_duration_*.R`                             |
| **Transition Models** | Define how states change                   | `R/models_transition_*.R`                           |
| **Factories**         | Build model components automatically       | `R/models_factory.R`                                |
| **Algorithms**        | EM iterations, log-likelihood, convergence | `R/algorithms_em.R`                                 |
| **C++ Core**          | Forward-backward & gamma computation       | `src/backward_forward.cpp`, `src/compute_gamma.cpp` |

---

## EM Algorithm

1. **E-step:**

   * Compute posterior probabilities via the forward–backward algorithm (`backward_forward.cpp`).
2. **M-step:**

   * Update emission, transition, and duration parameters using weighted likelihoods and GLMs.
3. **Convergence:**

   * Stop when the change in log-likelihood < `tol`.

---

## Example: HSMM with Covariates

```r
# Fit HSMM with wind speed as duration covariate
fit <- fit_hsmm(
  data = torus_data,
  n_states = 3,
  family = "toroidal",
  covariates_q = wind_speed,
  covariates_omega = temperature
)

# Extract dwell time distributions
fit$duration.model$get_params()
```

---


## 📚 References

* Lagona, F., & Mingione, M. (2025).
  *Nonhomogeneous Hidden Semi-Markov Models for Toroidal Data.*
  *JRSS Series C*, 74(1), 142–166. [DOI](https://doi.org/10.1093/jrsssc/qlae049)

* Mingione, M., Di Loro, P. A., Lagona, F., & Maruotti, A. (2025).
  *Environmental Risk Assessment via Nonhomogeneous Hidden Semi-Markov Models with Penalized VAR.*
  *arXiv preprint*, [arXiv:2509.14387v1](https://arxiv.org/abs/2509.14387v1)

---
