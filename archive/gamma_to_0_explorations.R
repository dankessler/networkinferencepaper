library('tidyverse')
library('latex2exp')
library('patchwork')
ggplot2::theme_set(theme_minimal())

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
h1_inv_deriv <- function(x, gamma) {
  c0 <- log(gamma / (1 - gamma))
  return((expit(logit(x) + c0) / (1 + exp(logit(x) + c0))) / (x*(1-x)))
}
h0_inv_deriv <- function(x, gamma) {
  c1 <- log((1 - gamma) / gamma)
  return((expit(logit(x) + c1) / (1 + exp(logit(x) + c1))) / (x*(1-x)))
}

# G2 piece functions (basically h1_inv_deriv and h0_inv_deriv without the 1 / (x*(1-x)))
G2_0_function <- function(x, gamma) {
  c0 <- log(gamma / (1 - gamma))
  return((expit(logit(x) - c0) / (1 + exp(logit(x) + c0))))
}
G2_1_function <- function(x, gamma) {
  c0 <- log(gamma / (1 - gamma))
  return((expit(logit(x) + c0) / (1 + exp(logit(x) - c0))))
}

# SIMULATION 1 - Histogram of variance as gamma gets small. We expect these
# to get very large since we think the variance should be blowing up. BTW this
# is the "TRUE" variance so nothing is being estimated here.

# ===========================
# == Simulation parameters ==
# ===========================

n <- 500                                        # Sample size
num_sim_per_gamma <- 300                  # Repetitions of the simulation
num_gamma_check <- 50
gamma_check <- seq(0.001, 0.50, length.out = num_gamma_check)
alpha <- 0.10

M <- runif(n, min = 0.05, max = 0.95)

# Where to store all the results of the simulation
D1_0_results <- array(numeric(), dim = c(num_gamma_check, num_sim_per_gamma))
G1_0_gamma_results <- array(numeric(), dim = c(num_gamma_check, num_sim_per_gamma))
G2_0_gamma_results <- array(numeric(), dim = c(num_gamma_check, num_sim_per_gamma))

D1_1_results <- array(numeric(), dim = c(num_gamma_check, num_sim_per_gamma))
G1_1_gamma_results <- array(numeric(), dim = c(num_gamma_check, num_sim_per_gamma))
G2_1_gamma_results <- array(numeric(), dim = c(num_gamma_check, num_sim_per_gamma))

variance_results <- array(numeric(), dim = c(num_gamma_check, num_sim_per_gamma))

for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]

  for (rep in 1:num_sim_per_gamma) {
    # Draw data and fission
    A <- rbinom(n, 1, M)
    Zfission <- rbinom(n, 1, gamma)
    Atr <- A * (1 - Zfission) + (1 - A) * Zfission

    # Indices for 0s and 1s in Atr
    A0_indices <- Atr == 0
    A1_indices <- !A0_indices
    num_0s <- sum(A0_indices)
    num_1s <- sum(A1_indices)
    num_all <- num_0s + num_1s

    # Conditional mean
    Cmask <- (gamma / (1-gamma))^(2*Atr - 1)
    Tvec <- M / (M + (1-M) * Cmask)

    B0_star <- mean(Tvec[A0_indices])
    B1_star <- mean(Tvec[A1_indices])

    # Variance calculations
    D1_0_results[gamma_index, rep] <- num_0s^2 / num_all^2
    D1_1_results[gamma_index, rep] <- num_1s^2 / num_all^2

    G1_0_gamma_results[gamma_index, rep] <- 1 / (num_0s * B0_star * (1 - B0_star))
    G1_1_gamma_results[gamma_index, rep] <- 1 / (num_1s * B1_star * (1 - B1_star))

    G2_0_gamma_results[gamma_index, rep] <- (G2_0_function(B0_star, gamma))^2
    G2_1_gamma_results[gamma_index, rep] <- (G2_1_function(B1_star, gamma))^2

    variance_results[gamma_index, rep] <-
      D1_0_results[gamma_index, rep] * G1_0_gamma_results[gamma_index, rep] * G2_0_gamma_results[gamma_index, rep] +
      D1_1_results[gamma_index, rep] * G1_1_gamma_results[gamma_index, rep] * G2_1_gamma_results[gamma_index, rep]
  }
}

# ---------------------
# -- Display results --
# ---------------------

gamma_display_indices <- 1:num_gamma_check

show_index <- 50

plot_df <- data.frame(variances = variance_results[gamma_display_indices[show_index], ])
ggplot(plot_df) +
  geom_histogram(aes(x = variances)) +
  ggtitle(paste0('Gamma = ', gamma_check[gamma_display_indices[show_index]]))

# Plot average variance by gamma
plot_df <- data.frame(gamma = gamma_check, avg_variance = apply(variance_results, 1, mean))
ggplot(plot_df) +
  geom_line(aes(x = gamma, y = avg_variance)) + xlab('Gamma') + ylab('Average variance')
