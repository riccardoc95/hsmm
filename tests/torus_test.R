library(R6)

# ---- SCRIPT DI TEST ----

set.seed(123)

# Parametri per test
n_obs <- 100
n_states <- 3
Pi <- 2 * pi  # valore per Pi che ti serve

# Genera dati random su (theta1, theta2) su [0, 2pi]
data <- cbind(runif(n_obs, 0, 2 * pi), runif(n_obs, 0, 2 * pi))

# Genera posteri probabilistici casuali (righe sommano a 1)
post_pi_raw <- matrix(runif(n_obs * n_states), nrow = n_obs, ncol = n_states)
post_pi <- post_pi_raw / rowSums(post_pi_raw)

# Prepara params
params <- list(
  data = data,
  n_states = n_states
)

# Crea istanza modello
torus_mod <- TorusModel$new(params, Pi, post_pi)

# Controlla parametri inizializzati
cat("Parametri iniziali:\n")
print(torus_mod$tau.mu1)
print(torus_mod$tau.mu2)
print(torus_mod$tau.k1)
print(torus_mod$tau.k2)
print(torus_mod$tau.rho)

# Controlla densità calcolata
cat("Densità iniziale dimensione:", dim(torus_mod$density), "\n")
print(head(torus_mod$density))

# Prova update con nuovi posteri casuali
post_pi_new_raw <- matrix(runif(n_obs * n_states), nrow = n_obs, ncol = n_states)
post_pi_new <- post_pi_new_raw / rowSums(post_pi_new_raw)

torus_mod$update(params, Pi, post_pi_new)

cat("Parametri dopo update:\n")
print(torus_mod$tau.mu1)
print(torus_mod$tau.mu2)
print(torus_mod$tau.k1)
print(torus_mod$tau.k2)
print(torus_mod$tau.rho)

cat("Densità dopo update dimensione:", dim(torus_mod$density), "\n")
print(head(torus_mod$density))

