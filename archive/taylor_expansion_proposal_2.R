library('networkinference')
library('latex2exp')
library('tidyverse')

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
h1 <- function(x, gamma) {
  return(expit(logit(x) + log((1 - gamma) / gamma)))
}
h0_inv <- function(x, gamma) {
  return(h1(x, gamma))
}
h1_inv <- function(x, gamma) {
  return(h0(x, gamma))
}
h0_deriv <- function(x, gamma) {
  r <- (1 - gamma) / gamma
  return((r / (x + (1-x)*r))^2)
}
h1_deriv <- function(x, gamma) {
  r <- gamma / (1 - gamma)
  return((r / (x + (1-x)*r))^2)
}

# ==========================================
# == See if a correction is even possible ==
# ==========================================

n <- 100
gamma_check <- seq(0.01, 0.50, length.out = 25)
num_sim_per_gamma <- 200
M <- runif(n, min = 0.1, max = 0.9)
alpha <- 0.10

for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]

  for (rep in 1:num_sim_per_gamma) {

    # Draw data and fission it
    A <- rbinom(n, size = 1, prob = M)
    W <- rbinom(n, size = 1, prob = gamma)
    A_tr <- A * (1 - W) + (1 - A) * W
    Tm <- A_tr * (h1(M, gamma)) + (1 - A_tr) * (h0(M, gamma))

    # Weights for zeros
    w0 <- sum(A_tr == 0) / n
    w1 <- 1 - w0

    # ------------------
    # -- Do inference --
    # ------------------

    # Varphi (what we actually care about)
    varphi <- mean(A)

    # Theta (what we don't care about but can target)
    theta <- mean(Tm)

    # Xi (Ethan's proposal)
    Phi0_star <- mean(Tm[as.logical(1-A_tr)])
    Phi1_star <- mean(Tm[as.logical(A_tr)])
    xi0 <- h0_inv(Phi0_star, gamma)
    xi1 <- h1_inv(Phi1_star, gamma)
    xi <- w0 * xi0 + w1 * xi1

    # Capricorn (Anna's proposal)
    ci_offsets <- (2*A_tr - 1) * log((1-gamma) / gamma)
    capricorn <- expit(glm(A ~ 1 + offset(ci_offsets), family = 'binomial')$coefficients)

    # Aquarius (Daniela's Taylor expansion proposal)


    # ====================
    # == Record coverge ==
    # ====================

    varphi_record[gamma_index, rep] <- varphi

    theta_record[gamma_index, rep] <- theta
    theta_hat_record[gamma_index, rep] <- theta_hat
    theta_ci_lower_record[gamma_index, rep] <- theta_hat - theta_moe
    theta_ci_upper_record[gamma_index, rep] <- theta_hat + theta_moe
    theta_covered_record[gamma_index, rep] <- (theta_hat - theta_moe <= theta) & (theta <= theta_hat + theta_moe)

    xi_record[gamma_index, rep] <- xi
    xi_hat_record[gamma_index, rep] <- xi_hat
    xi_ci_lower_record[gamma_index, rep] <- xi_hat - xi_moe
    xi_ci_upper_record[gamma_index, rep] <- xi_hat + xi_moe
    xi_covered_record[gamma_index, rep] <- (xi_hat - xi_moe <= xi) & (xi <= xi_hat + xi_moe)

    varphi_covered_by_xi_record[gamma_index, rep] <- (xi_hat - xi_moe <= varphi) & (varphi <= xi_hat + xi_moe)
    varphi_covered_by_theta_record[gamma_index, rep] <- (theta_hat - theta_moe <= varphi) & (varphi <= theta_hat + theta_moe)

    Capricorn_record[gamma_index, rep] <- Capricorn
  }
}
