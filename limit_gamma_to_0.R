library(tidyverse)

# =======================
# == Helpful functions ==
# =======================

logit <- function(x) {
  return(log(x / (1 - x)))
}
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}
h0 <- function(x, gamma) {
  return(expit(logit(x) + log(gamma / (1 - gamma))))
}
h0_alt <- function(x, gamma) {
  return(x / (x + (1-x)*((1-gamma) / gamma)))
}
h1 <- function(x, gamma) {
  return(expit(logit(x) + log((1 - gamma) / gamma)))
}
h0_inv <- function(x, gamma) {
  return(h1(x, gamma))
}
h1_inv <- function(x, gamma) {
  return(h0(x, gamma))
}
# h1_inv_deriv <- function(x, gamma) {
#   c0 <- log(gamma / (1 - gamma))
#   return((expit(logit(x) + c0) / (1 + exp(logit(x) + c0))) / (x*(1-x)))
# }
h1_inv_deriv <- function(x, gamma) {
  c0 <- log(gamma / (1 - gamma))
  return((expit(logit(x) + c0) / (1 + exp(logit(x) + c0))) / (x*(1-x)))
}
h0_inv_deriv <- function(x, gamma) {
  c1 <- log((1 - gamma) / gamma)
  return((expit(logit(x) + c1) / (1 + exp(logit(x) + c1))) / (x*(1-x)))
}

# Potential limiting things
# NEEDS TO BE UPDATED TO ACCOUNT FORGOTTEN 1/(X*(1-X)) FACTOR
A2_0_limit <- function(m, gamma) {
  return((mean((m / (1-m)))^2) / ((1 + mean(m / (1-m)))^4))
}


# SIM 1
# ------

# Simulation parameters
set.seed(1)
n <- 300
n_rep <- 150

M <- runif(n, min = 0.01, max = 0.99)
gamma_check <- seq(0.001, 0.50, length.out = 100)

# Stored results
variance_estimates <- array(numeric(), c(length(gamma_check), n_rep))
variance_limits_suspected <- array(numeric(), c(length(gamma_check), n_rep))

# Stored results 2
a2_actual <- rep(NA, n_rep)
a2_limiting <- rep(NA, n_rep)

for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]

  for (rep in 1:n_rep) {
    # Draw data and fission
    A <- rbinom(n, 1, M)
    Zfission <- rbinom(n, 1, gamma)
    Atr <- A * (1 - Zfission) + (1 - A) * Zfission

    # Indices for 0s and 1s in Atr
    A0_indices <- Atr == 0
    A1_indices <- !A0_indices

    # Conditional mean
    Cmask <- (gamma / (1-gamma))^(2*Atr - 1)
    Tvec <- M / (M + (1-M) * Cmask)

    Lambda0_star <- mean(Tvec[A0_indices])
    Lambda1_star <- mean(Tvec[A1_indices])

    # Variance estimate
    var_est <- h0_inv_deriv(Lambda0_star, gamma)^2

    # Suspected variance limit
    var_limit_suspect <- A2_0_limit(M[A0_indices], gamma)

    # Store everything
    variance_estimates[gamma_index, rep] <- var_est
    variance_limits_suspected[gamma_index, rep] <- var_limit_suspect
  }
}

var_est_avg <- apply(variance_estimates, 1, mean)
var_limit_avg <- apply(variance_limits_suspected, 1, mean)

# Plot
plot_df <- data.frame(gamma = gamma_check,
                      var_est_avg = var_est_avg,
                      var_limit_avg = var_limit_avg)
ggplot(plot_df) +
  geom_line(aes(x = gamma, y = var_est_avg)) +
  geom_line(aes(x = gamma, y = var_limit_avg))










# SIM 2

# Simulation parameters
set.seed(1)
n <- 200
n_rep <- 2000

M <- runif(n, min = 0.01, max = 0.99)
gamma <- 0.0001

# Stored results 2
a2_actual <- rep(NA, n_rep)
a2_limiting <- rep(NA, n_rep)

for (rep in 1:n_rep) {
  # Draw data and fission
  A <- rbinom(n, 1, M)
  Zfission <- rbinom(n, 1, gamma)
  Atr <- A * (1 - Zfission) + (1 - A) * Zfission

  # Indices for 0s and 1s in Atr
  A0_indices <- Atr == 0
  A1_indices <- !A0_indices

  # Conditional mean
  Cmask <- (gamma / (1-gamma))^(2*Atr - 1)
  Tvec <- M / (M + (1-M) * Cmask)

  Lambda0_star <- mean(Tvec[A0_indices])
  Lambda1_star <- mean(Tvec[A1_indices])

  # Variance estimate
  var_est <- h0_inv_deriv(Lambda0_star, gamma)^2

  # Suspected variance limit
  var_limit_suspect <- A2_0_limit(M[A0_indices], gamma)

  # Store everything
  a2_actual[rep] <- var_est
  a2_limiting[rep] <- var_limit_suspect
}

# Plot
plot_df <- data.frame(a2_actual = a2_actual,
                      a2_limiting = a2_limiting)
ggplot(plot_df) +
  geom_point(aes(x = a2_actual, y = a2_limiting))

