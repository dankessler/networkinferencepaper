# ============================================
# == Simulations for network thinning paper ==
# == (Binary edges)                         ==
# == Author: Ethan Ancell                   ==
# ============================================

# Libraries
library(ggplot2)
library(tidyverse)
library(nett)
library(mclust)
library(latex2exp)

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

# =========================
# == Simulation settings ==
# =========================

set.seed(1)

# Most comprehensive settings
# (Currently commented out)
# ---------------------------

# n_check <- c(50, 100, 200)                   # Total network size
# n_per_community_check <- c(1, 2, 10)         # Amount of nodes in each true community
# K_check <- c(2, 5, 10)                       # Guessed number of communities
# num_sim <- 500                               # Number of replications
# signal_regimes <- list(c(11, 10), c(15, 10)) # Signal regimes for mean matrix
# gamma_check <- c(0.2, 0.5, 0.7)


# Feel free to adjust any of these
# settings on the fly for whatever you're doing
# ---------------------------------------------

# -------------
# -- Table 1 --
# -------------

set.seed(1)
num_sim <- 5000 # Number of replications

n_check <- c(100, 200, 500) # Total network size
K_check <- c(2, 5, 10) # Guessed number of communities
K_true_check <- c(5)
gamma_check <- c(0.25)
signal_regimes <- list(c(0.75, 0.5)) # Signal regimes for mean matrix

alpha <- 0.10
use_random_u <- FALSE

# --------------
# -- Figure 1 --
# --------------

set.seed(1)
num_sim <- 5000 # Number of replications

n_check <- c(200) # Total network size
K_check <- c(2, 5, 10) # Guessed number of communities
K_true_check <- c(5)
gamma_check <- c(0.25)
signal_regimes <- list(c(0.75, 0.5)) # Signal regimes for mean matrix

alpha <- 0.10
use_random_u <- FALSE

# --------------
# -- Figure 2 --
# --------------

set.seed(1)
num_sim <- 5000 # Number of replications

n_check <- c(200) # Total network size
K_check <- c(5) # Guessed number of communities
K_true_check <- c(5)
gamma_check <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)
signal_regimes <- list(c(0.75, 0.55), c(0.75, 0.50), c(0.75, 0.45), c(0.75, 0.40), c(0.75, 0.35)) # Signal regimes for mean matrix

alpha <- 0.10
use_random_u <- FALSE

# =======================
# == Set up simulation ==
# =======================

# Useful for bookkeeping and plotting: which dimensions correspond to which
# things that we are looping through.
looping_dimensions <- c('n', 'K_true', 'K', 'Signal Regime')

# Unpack rho_1 and rho_2 from the signal regimes
rho_1 <- signal_regimes[[1]][1] # Assumes that this stays the same across all of them
rho_2_check <- rep(NA, length(signal_regimes))
for (i in 1:length(signal_regimes)) {
  rho_2_check[i] <- signal_regimes[[i]][2]
}
rho_1_minus_rho_2 <- rho_1 - rho_2_check

# Where results are stored
# ------------------------

# Estimates, targets, variances, etc.
# (saving all of this makes it easier to calculate some things without
# re-running the entire experiment)
xi_est_thinning <- array(0, dim = c(length(n_check),
                                    length(K_true_check),
                                    length(K_check),
                                    length(signal_regimes),
                                    length(gamma_check),
                                    num_sim))
xi_target_thinning <- array(0, dim = c(length(n_check),
                                       length(K_true_check),
                                       length(K_check),
                                       length(signal_regimes),
                                       length(gamma_check),
                                       num_sim))
theta_target_thinning <- array(0, dim = c(length(n_check),
                                          length(K_true_check),
                                          length(K_check),
                                          length(signal_regimes),
                                          length(gamma_check),
                                          num_sim))
xi_variances_thinning <- array(0, dim = c(length(n_check),
                                       length(K_true_check),
                                       length(K_check),
                                       length(signal_regimes),
                                       length(gamma_check),
                                       num_sim))

theta_est_naive <- array(0, dim = c(length(n_check),
                                    length(K_true_check),
                                    length(K_check),
                                    length(signal_regimes),
                                    length(gamma_check),
                                    num_sim))
theta_target_naive <- array(0, dim = c(length(n_check),
                                  length(K_true_check),
                                  length(K_check),
                                  length(signal_regimes),
                                  length(gamma_check),
                                  num_sim))
theta_variances_naive <- array(0, dim = c(length(n_check),
                                    length(K_true_check),
                                    length(K_check),
                                    length(signal_regimes),
                                    length(gamma_check),
                                    num_sim))
rand_results <- array(0, dim = c(length(n_check),
                                 length(K_true_check),
                                 length(K_check),
                                 length(signal_regimes),
                                 length(gamma_check),
                                 num_sim))

# Coverage results of (1-alpha) confidence intervals
xi_covered_thinning <- array(0, dim = c(length(n_check),
                                        length(K_true_check),
                                        length(K_check),
                                        length(signal_regimes),
                                        length(gamma_check)))
theta_covered_thinning <- array(0, dim = c(length(n_check),
                                           length(K_true_check),
                                           length(K_check),
                                           length(signal_regimes),
                                           length(gamma_check)))
theta_covered_naive <- array(0, dim = c(length(n_check),
                                  length(K_true_check),
                                  length(K_check),
                                  length(signal_regimes),
                                  length(gamma_check)))

# ==============
# == Simulate ==
# ==============

# Sample size
for (n_index in 1:length(n_check)) {
  n <- n_check[n_index]
  cat(paste0('(n = ', n, ')\n'))

  # True communities
  for (K_true_index in 1:length(K_true_check)) {
    # n_per_community <- n_per_community_check[n_per_community_index]
    K_true <- K_true_check[K_true_index]
    n_per_community <- n / K_true
    cat(paste0('\t(K_true = ', K_true, ')\n'))

    # Number of communities used in clustering
    for (K_index in 1:length(K_check)) {
      K <- K_check[K_index]
      cat(paste0('\t\t(K = ', K, ')\n'))

      # Signal regimes in mean matrix
      for (signal_regime_index in 1:length(signal_regimes)) {
        rho1 <- signal_regimes[[signal_regime_index]][1]
        rho2 <- signal_regimes[[signal_regime_index]][2]
        cat(paste0('\t\t\t(rho1 = ', rho1, ', rho2 = ', rho2, ')\n'))

        # Different values of gamma in fission
        for (gamma_index in 1:length(gamma_check)) {
          gamma <- gamma_check[gamma_index]
          cat(paste0('\t\t\t\t(gamma = ', gamma, ')\n'))

          # -------------------------------------------------
          # -- Actual simulation part (non just for-loops) --
          # -------------------------------------------------

          # Create mean matrices
          Theta_true <- matrix(rep(rho2, K_true^2), nrow = K_true)
          diag(Theta_true) <- rho1
          M <- matrix(rep(NA, n^2), nrow = n)
          for (k in 1:K_true) {
            for (l in 1:K_true) {
              M[((k-1)*n_per_community + 1):(k*n_per_community),
                ((l-1)*n_per_community + 1):(l*n_per_community)] <- Theta_true[k, l]
            }
          }
          # "True" community membership
          z_true <- rep(1:K_true, each = n_per_community)

          for (rep in 1:num_sim) {

            # Sample u
            if (use_random_u) {
              # Sample the vector u from the unit K^2-sphere
              u <- rnorm(K^2)
              u <- u / sqrt(sum(u^2))
            } else {
              u <- c(1, rep(0, K^2 - 1))
            }


            # The purpose of maybe needing to redraw a matrix is that sometimes
            # you get a single node in an estimated community, and we want to
            # disallow that just for the sake of simulation, because then you
            # get infinite variance estimates and it gets a little nasty.

            need_good_matrix <- TRUE
            while (need_good_matrix) {

              # Draw A
              A <- matrix(rbinom(n = n^2, size = 1, prob = as.vector(M)), nrow = n)

              # Fission A
              Zfission <- matrix(rbinom(n = n^2, size = 1, prob = gamma), nrow = n)
              A_tr <- A * (1 - Zfission) + (1 - A) * Zfission

              # Cluster (thinning)
              z_hat_initial <- nett::spec_clust(A_tr, K = K)
              z_hat <- matrix(rep(NA, n*K), nrow = n)
              for (i in 1:K) {
                z_hat[, i] <- 1 * (z_hat_initial == i)
              }
              n_hat <- apply(z_hat, 2, sum)

              # Check RAND index agreement between true clusters and non true clusters
              rand_results[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
                mclust::adjustedRandIndex(z_hat_initial, z_true)

              # Cluster (naive)
              # z_hat_initial_naive <- nett::spec_clust(A_2, K = K)
              z_hat_initial_naive <- nett::spec_clust(A, K = K)
              z_hat_naive <- matrix(rep(NA, n*K), nrow = n)
              for (i in 1:K) {
                z_hat_naive[, i] <- 1 * (z_hat_initial_naive == i)
              }
              n_hat_naive <- apply(z_hat_naive, 2, sum)

              # ---------------------
              # -- Useful matrices --
              # ---------------------

              NN_inv <- diag(1 / diag(t(z_hat) %*% z_hat))
              NN_inv_naive <- diag(1 / diag(t(z_hat_naive) %*% z_hat_naive))

              comm_pair_sample_size <- t(z_hat) %*% matrix(1, nrow = n, ncol = n) %*% z_hat
              comm_pair_sample_size_ones <- t(z_hat) %*% A_tr %*% z_hat
              comm_pair_sample_size_zeros <- t(z_hat) %*% (1 - A_tr) %*% z_hat

              comm_pair_sample_size_naive <- t(z_hat_naive) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_naive

              # Weighting of ones and zeros
              ones_weighting <- comm_pair_sample_size_ones / comm_pair_sample_size
              zeros_weighting <- comm_pair_sample_size_zeros / comm_pair_sample_size

              # --------------------------------------
              # -- Targets/estimates/variances/etc. --
              # --------------------------------------

              # Thinning
              # --------

              Cmask <- (gamma / (1-gamma))^(2*A_tr - 1) # The constant gamma/(1-gamma) masking thing
              Tmat <- M / (M + (1-M)*Cmask) # Calculate T as a function of M

              # Estimators and estimands related to xi(Atr)
              Lambda0_star <- (t(z_hat) %*% (Tmat * (1 - A_tr)) %*% z_hat) /
                comm_pair_sample_size_zeros
              Lambda1_star <- (t(z_hat) %*% (Tmat * A_tr) %*% z_hat) /
                comm_pair_sample_size_ones
              Lambda0_hat <- (t(z_hat) %*% (A * (1 - A_tr)) %*% z_hat) /
                comm_pair_sample_size_zeros
              Lambda1_hat <- (t(z_hat) %*% (A * A_tr) %*% z_hat) /
                comm_pair_sample_size_ones

              Lambda0_var_ideal <- (t(z_hat) %*% (Tmat * (1 - Tmat) * (1 - A_tr)) %*% z_hat) /
                comm_pair_sample_size_zeros^2
              Lambda0_var_est <- (Lambda0_hat * (1 - Lambda0_hat)) / comm_pair_sample_size_zeros

              Lambda1_var_ideal <- (t(z_hat) %*% (Tmat * (1 - Tmat) * A_tr) %*% z_hat) /
                comm_pair_sample_size_ones^2
              Lambda1_var_est <- (Lambda1_hat * (1 - Lambda1_hat)) / comm_pair_sample_size_ones

              Phi_matrix <- zeros_weighting * h0_inv(Lambda0_star, gamma) + ones_weighting * h1_inv(Lambda1_star, gamma)
              Phi_hat_matrix <- zeros_weighting * h0_inv(Lambda0_hat, gamma) + ones_weighting * h1_inv(Lambda1_hat, gamma)

              Phi_var_matrix_ideal <- (zeros_weighting^2 * (Lambda0_var_ideal * (h0_inv_deriv(Lambda0_star, gamma))^2)) +
                (ones_weighting^2 * (Lambda1_var_ideal * (h1_inv_deriv(Lambda1_star, gamma))^2))

              Phi_var_matrix_est <- zeros_weighting^2 * Lambda0_var_est * (h0_inv_deriv(Lambda0_hat, gamma))^2 +
                ones_weighting^2 * Lambda1_var_est * (h1_inv_deriv(Lambda1_hat, gamma))^2

              xi <- t(u) %*% as.vector(Phi_matrix)
              theta <- t(u) %*% as.vector(NN_inv %*% t(z_hat) %*% M %*% z_hat %*% NN_inv)
              xi_hat <- t(u) %*% as.vector(Phi_hat_matrix)
              xi_moe_ideal <- qnorm(1 - alpha/2) * sqrt(t(u) %*% diag(as.vector(Phi_var_matrix_ideal)) %*% u)
              xi_moe_est <- qnorm(1 - alpha/2) * sqrt(t(u) %*% diag(as.vector(Phi_var_matrix_est)) %*% u)

              # Naive
              # -----
              theta_naive_matrix <- NN_inv_naive %*% t(z_hat_naive) %*% M %*% z_hat_naive %*% NN_inv_naive
              theta_est_naive_matrix <- NN_inv_naive %*% t(z_hat_naive) %*% A %*% z_hat_naive %*% NN_inv_naive

              theta_naive <- t(u) %*% as.vector(theta_naive_matrix)
              theta_est_naive_temp <- t(u) %*% as.vector(theta_est_naive_matrix)

              theta_var_naive_matrix <- (theta_est_naive_matrix * (1 - theta_est_naive_matrix)) / comm_pair_sample_size_naive
              theta_var_naive <- t(u) %*% diag(as.vector(theta_var_naive_matrix)) %*% u
              theta_moe_naive <- qnorm(1 - alpha / 2) * sqrt(theta_var_naive)

              # We're good to go if there are no community pairs with just one
              # member in it in the Atr = 0 and Atr = 1 splits, and also if
              # there are no NaNs which appear in the estimates (which happens if
              # there are all 1s or all 0s inside of a community pair.)
              if ((sum(comm_pair_sample_size_ones == 1) == 0) &
                  (sum(comm_pair_sample_size_zeros == 1) == 0) &
                  (sum(comm_pair_sample_size_naive == 1) == 0) &
                  (!is.nan(xi_moe_est)) &
                  (!is.nan(theta_moe_naive))) {
                need_good_matrix <- FALSE
              }
            }

            # Save results
            # ------------
            xi_est_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              xi_hat
            xi_target_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              xi
            theta_target_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              theta
            xi_variances_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              t(u) %*% diag(as.vector(Phi_var_matrix_est)) %*% u

            theta_est_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              theta_est_naive_temp
            theta_target_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              theta_naive
            theta_variances_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              theta_var_naive

            # Coverage results of (1-alpha) confidence intervals
            xi_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] <-
              xi_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] +
              as.integer((xi_hat - xi_moe_est <= xi) & (xi <= xi_hat + xi_moe_est))
            theta_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] <-
              theta_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] +
              as.integer((xi_hat - xi_moe_est <= theta) & (theta <= xi_hat + xi_moe_est))

            theta_covered_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index] <-
              theta_covered_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index] +
              as.integer((theta_est_naive_temp - theta_moe_naive <= theta_naive) & (theta_naive <= theta_est_naive_temp + theta_moe_naive))
          }
        }
      }
    }
  }
}

xi_covered_thinning <- xi_covered_thinning / num_sim
theta_covered_thinning <- theta_covered_thinning / num_sim
theta_covered_naive <- theta_covered_naive / num_sim

xi_covered_thinning
theta_covered_thinning
theta_covered_naive

# ====================================================
# == SAVE the simulation results after running them ==
# == Uncomment/comment the lines appropriately.     ==
# ====================================================

# -----------
# - Table 1 -
# -----------

# saveRDS(xi_est_thinning, file = 'save_data/SIM1_BERNOULLI_xi_est_thinning.RDS')
# saveRDS(xi_target_thinning, file = 'save_data/SIM1_BERNOULLI_xi_target_thinning.RDS')
# saveRDS(theta_target_thinning, file = 'save_data/SIM1_BERNOULLI_theta_target_thinning.RDS')
# saveRDS(xi_variances_thinning, file = 'save_data/SIM1_BERNOULLI_xi_variances_thinning.RDS')
# saveRDS(theta_est_naive, file = 'save_data/SIM1_BERNOULLI_theta_est_naive.RDS')
# saveRDS(theta_target_naive, file = 'save_data/SIM1_BERNOULLI_theta_target_naive.RDS')
# saveRDS(theta_variances_naive, file = 'save_data/SIM1_BERNOULLI_theta_variances_naive.RDS')
# saveRDS(rand_results, file = 'save_data/SIM1_BERNOULLI_rand_results.RDS')
# saveRDS(xi_covered_thinning, file = 'save_data/SIM1_BERNOULLI_xi_covered_thinning.RDS')
# saveRDS(theta_covered_thinning, file = 'save_data/SIM1_BERNOULLI_theta_covered_thinning.RDS')
# saveRDS(theta_covered_naive, file = 'save_data/SIM1_BERNOULLI_theta_covered_naive.RDS')

xi_est_thinning <- readRDS('save_data/SIM1_BERNOULLI_xi_est_thinning.RDS')
xi_target_thinning <- readRDS('save_data/SIM1_BERNOULLI_xi_target_thinning.RDS')
theta_target_thinning <- readRDS('save_data/SIM1_BERNOULLI_theta_target_thinning.RDS')
xi_variances_thinning <- readRDS('save_data/SIM1_BERNOULLI_xi_variances_thinning.RDS')
theta_est_naive <- readRDS('save_data/SIM1_BERNOULLI_theta_est_naive.RDS')
theta_target_naive <- readRDS('save_data/SIM1_BERNOULLI_theta_target_naive.RDS')
theta_variances_naive <- readRDS('save_data/SIM1_BERNOULLI_theta_variances_naive.RDS')
rand_results <- readRDS('save_data/SIM1_BERNOULLI_rand_results.RDS')
xi_covered_thinning <- readRDS('save_data/SIM1_BERNOULLI_xi_covered_thinning.RDS')
theta_covered_thinning <- readRDS('save_data/SIM1_BERNOULLI_theta_covered_thinning.RDS')
theta_covered_naive <- readRDS('save_data/SIM1_BERNOULLI_theta_covered_naive.RDS')

# ------------
# - Figure 1 -
# ------------

# saveRDS(xi_est_thinning, file = 'save_data/SIM2_BERNOULLI_xi_est_thinning.RDS')
# saveRDS(xi_target_thinning, file = 'save_data/SIM2_BERNOULLI_xi_target_thinning.RDS')
# saveRDS(theta_target_thinning, file = 'save_data/SIM2_BERNOULLI_theta_target_thinning.RDS')
# saveRDS(xi_variances_thinning, file = 'save_data/SIM2_BERNOULLI_xi_variances_thinning.RDS')
# saveRDS(theta_est_naive, file = 'save_data/SIM2_BERNOULLI_theta_est_naive.RDS')
# saveRDS(theta_target_naive, file = 'save_data/SIM2_BERNOULLI_theta_target_naive.RDS')
# saveRDS(theta_variances_naive, file = 'save_data/SIM2_BERNOULLI_theta_variances_naive.RDS')
# saveRDS(rand_results, file = 'save_data/SIM2_BERNOULLI_rand_results.RDS')
# saveRDS(xi_covered_thinning, file = 'save_data/SIM2_BERNOULLI_xi_covered_thinning.RDS')
# saveRDS(theta_covered_thinning, file = 'save_data/SIM2_BERNOULLI_theta_covered_thinning.RDS')
# saveRDS(theta_covered_naive, file = 'save_data/SIM2_BERNOULLI_theta_covered_naive.RDS')

xi_est_thinning <- readRDS('save_data/SIM2_BERNOULLI_xi_est_thinning.RDS')
xi_target_thinning <- readRDS('save_data/SIM2_BERNOULLI_xi_target_thinning.RDS')
theta_target_thinning <- readRDS('save_data/SIM2_BERNOULLI_theta_target_thinning.RDS')
xi_variances_thinning <- readRDS('save_data/SIM2_BERNOULLI_xi_variances_thinning.RDS')
theta_est_naive <- readRDS('save_data/SIM2_BERNOULLI_theta_est_naive.RDS')
theta_target_naive <- readRDS('save_data/SIM2_BERNOULLI_theta_target_naive.RDS')
theta_variances_naive <- readRDS('save_data/SIM2_BERNOULLI_theta_variances_naive.RDS')
rand_results <- readRDS('save_data/SIM2_BERNOULLI_rand_results.RDS')
xi_covered_thinning <- readRDS('save_data/SIM2_BERNOULLI_xi_covered_thinning.RDS')
theta_covered_thinning <- readRDS('save_data/SIM2_BERNOULLI_theta_covered_thinning.RDS')
theta_covered_naive <- readRDS('save_data/SIM2_BERNOULLI_theta_covered_naive.RDS')

# ------------
# - Figure 2 -
# ------------

# saveRDS(xi_est_thinning, file = 'save_data/SIM3_BERNOULLI_xi_est_thinning.RDS')
# saveRDS(xi_target_thinning, file = 'save_data/SIM3_BERNOULLI_xi_target_thinning.RDS')
# saveRDS(theta_target_thinning, file = 'save_data/SIM3_BERNOULLI_theta_target_thinning.RDS')
# saveRDS(xi_variances_thinning, file = 'save_data/SIM3_BERNOULLI_xi_variances_thinning.RDS')
# saveRDS(theta_est_naive, file = 'save_data/SIM3_BERNOULLI_theta_est_naive.RDS')
# saveRDS(theta_target_naive, file = 'save_data/SIM3_BERNOULLI_theta_target_naive.RDS')
# saveRDS(theta_variances_naive, file = 'save_data/SIM3_BERNOULLI_theta_variances_naive.RDS')
# saveRDS(rand_results, file = 'save_data/SIM3_BERNOULLI_rand_results.RDS')
# saveRDS(xi_covered_thinning, file = 'save_data/SIM3_BERNOULLI_xi_covered_thinning.RDS')
# saveRDS(theta_covered_thinning, file = 'save_data/SIM3_BERNOULLI_theta_covered_thinning.RDS')
# saveRDS(theta_covered_naive, file = 'save_data/SIM3_BERNOULLI_theta_covered_naive.RDS')

xi_est_thinning <- readRDS('save_data/SIM3_BERNOULLI_xi_est_thinning.RDS')
xi_target_thinning <- readRDS('save_data/SIM3_BERNOULLI_xi_target_thinning.RDS')
theta_target_thinning <- readRDS('save_data/SIM3_BERNOULLI_theta_target_thinning.RDS')
xi_variances_thinning <- readRDS('save_data/SIM3_BERNOULLI_xi_variances_thinning.RDS')
theta_est_naive <- readRDS('save_data/SIM3_BERNOULLI_theta_est_naive.RDS')
theta_target_naive <- readRDS('save_data/SIM3_BERNOULLI_theta_target_naive.RDS')
theta_variances_naive <- readRDS('save_data/SIM3_BERNOULLI_theta_variances_naive.RDS')
rand_results <- readRDS('save_data/SIM3_BERNOULLI_rand_results.RDS')
xi_covered_thinning <- readRDS('save_data/SIM3_BERNOULLI_xi_covered_thinning.RDS')
theta_covered_thinning <- readRDS('save_data/SIM3_BERNOULLI_theta_covered_thinning.RDS')
theta_covered_naive <- readRDS('save_data/SIM3_BERNOULLI_theta_covered_naive.RDS')

# ==================================
# == Plotting and viewing results ==
# ==================================

# Table 1 - Checking coverage of thinning and naive
# -------------------------------------------------
plot_table_1 <- function(data) {

  signal_index <- 1
  K_true_index <- 1
  gamma_index <- 1

  # Print header column
  cat("|         |")
  for (K_index in 1:length(K_check)) {
    K <- K_check[K_index]
    cat(sprintf(" K = %03d |", K))
  }
  cat('\n')

  # Print rows of table below header
  for (n_index in 1:length(n_check)) {
    n <- n_check[n_index]

    cat(sprintf("| n = %03d |", n))

    for (K_index in 1:length(K_check)) {
      K <- K_check[K_index]

      cat(sprintf("   %.2f |", 100 * data[n_index, K_true_index, K_index, signal_index, gamma_index]))
    }
    cat('\n')
  }
}

cat('Naive\n')
plot_table_1(theta_covered_naive)

cat('Xi - Thinning\n')
plot_table_1(xi_covered_thinning)

cat('Theta - Thinning\n')
plot_table_1(theta_covered_thinning)

# Figure 1 - Showing a comparison of nominal vs empirical coverage
# ----------------------------------------------------------------

checked_coverages <- seq(0, 1.0, length.out = 100)
n_coverages <- length(checked_coverages)

n_index <- 1
K_true_index <- 1
K_true <- K_true_check[K_true_index]
signal_index <- 1

# Results
figure1_coverage_thinning <- array(numeric(), c(length(K_check), length(gamma_check), length(checked_coverages)))
figure1_coverage_naive <- array(numeric(), c(length(K_check), length(gamma_check), length(checked_coverages)))

for (K_index in 1:length(K_check)) {

  for (gamma_index in 1:length(gamma_check)) {

    for (coverage_index in 1:length(checked_coverages)) {

      coverage <- checked_coverages[coverage_index]
      z_quantile <- qnorm(1 - (1 - coverage) / 2)

      # Extract relevant slices from the results array
      targets_thinning_slice <- theta_target_thinning[n_index, K_true_index, K_index, signal_index, gamma_index, ]
      theta_target_naive_slice <- theta_target_naive[n_index, K_true_index, K_index, signal_index, gamma_index, ]

      # NOTE: This is currently set to do coverage for theta(Atr) so that it's "comparable" to theta(A)
      # in the naive approach.

      upper_bounds_thinning <-
        xi_est_thinning[n_index, K_true_index, K_index, signal_index, gamma_index, ] +
        z_quantile * sqrt(xi_variances_thinning[n_index, K_true_index, K_index, signal_index, gamma_index, ])
      lower_bounds_thinning <-
        xi_est_thinning[n_index, K_true_index, K_index, signal_index, gamma_index, ] -
        z_quantile * sqrt(xi_variances_thinning[n_index, K_true_index, K_index, signal_index, gamma_index, ])

      upper_bounds_naive <-
        theta_est_naive[n_index, K_true_index, K_index, signal_index, gamma_index, ] +
        z_quantile * sqrt(theta_variances_naive[n_index, K_true_index, K_index, signal_index, gamma_index, ])
      lower_bounds_naive <-
        theta_est_naive[n_index, K_true_index, K_index, signal_index, gamma_index, ] -
        z_quantile * sqrt(theta_variances_naive[n_index, K_true_index, K_index, signal_index, gamma_index, ])

      # Store coverage
      figure1_coverage_thinning[K_index, gamma_index, coverage_index] <-
        mean((lower_bounds_thinning <= targets_thinning_slice) & (targets_thinning_slice <= upper_bounds_thinning))
      figure1_coverage_naive[K_index, gamma_index, coverage_index] <-
        mean((lower_bounds_naive <= theta_target_naive_slice) & (theta_target_naive_slice <= upper_bounds_naive))

    }

  }

}

# Stack all this information in a data frame that can be used for plotting
plot_df <- data.frame()
for (K_index in 1:length(K_check)) {
  for (gamma_index in 1:length(gamma_check)) {

    # Create a method name for thinning / naive that includes K = ---
    # thinning_method_name <- paste0('Thinning (K = ', K_check[K_index], ')')
    # naive_method_name <- paste0('Naive (K = ', K_check[K_index], ')')

    # Thinning
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_thinning[K_index, gamma_index, ],
                                method = rep('Proposed', n_coverages),
                                # method_print = rep(thinning_method_name, n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                gamma = rep(gamma_check[gamma_index], n_coverages)))

    # Naive
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_naive[K_index, gamma_index, ],
                                method = rep('Naive', n_coverages),
                                # method_print = rep(naive_method_name, n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                gamma = rep(gamma_check[gamma_index], n_coverages)))
  }
}

# Create plot
legend_colors <- c('2' = 'cadetblue3', '5' = 'royalblue3', '10' = 'blue4')
figure1 <- ggplot()
for (K_index in 1:length(K_check)) {
  K <- K_check[K_index]

  # Subset data
  df_subset_thin <- plot_df %>% filter(method == 'Proposed', K == K)
  df_subset_naive <- plot_df %>% filter(method == 'Naive', K == K)

  # Add it onto the plot
  figure1 <- figure1 +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Proposed'), size = 0.9, data = df_subset_thin)
  figure1 <- figure1 +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Naive'), size = 0.7, data = df_subset_naive)
}
figure1 <- figure1 +
  scale_color_manual(breaks = c('2', '5', '10'),
                     values = legend_colors) +
  scale_linetype_manual(breaks = c('Proposed', 'Naive'),
                        values = c('Proposed' = 'solid', 'Naive' =
                                     'dashed')) +
  xlab('Nominal Coverage') + ylab('Empirical Coverage') +
  coord_fixed() +
  labs(color = 'K', linetype = "Method") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

figure1
ggsave('figures/nominal_empirical_coverage_bernoulli.pdf',
       plot = figure1, device = 'pdf', width = 4.5, height = 3.5)

# Figure 2 - Tradeoff of signal regime, gamma, RAND index, CI width
# -------------------------------------------------------------------

n_index <- 1
K_true_index <- 1
K_index <- 1
alpha <- 0.10

# Figure 2a - CI width as a function of signal, separate by gamma
# -----------------------------------------------------------------
ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(xi_variances_thinning[n_index, K_true_index, K_index, , , ])
# ci_width <- (2 * qnorm(1 - alpha / 2) * sqrt(xi_variances_thinning[n_index, K_true_index, K_index, , , ])) / sqrt(xi_est_thinning[n_index, K_true_index, K_index, , , ])
ci_width <- apply(ci_width, c(1, 2), mean)
gamma_index_plotting <- matrix(rep(gamma_check, length(signal_regimes)), ncol = length(gamma_check), byrow = TRUE)
rho_2_plotting <- matrix(rep(rho_2_check, length(gamma_check)), ncol = length(gamma_check))
rho_1_minus_rho_2_plotting <- matrix(rep(rho_1_minus_rho_2, length(gamma_check)), ncol = length(gamma_check))
signal_index_plotting <- matrix(rep(1:length(signal_regimes), length(gamma_check)), ncol = length(gamma_check))

plot_df <- data.frame(signal_index = as.vector(signal_index_plotting),
                      rho_2 = as.vector(rho_2_plotting),
                      rho_1_minus_rho_2 = as.character(as.vector(rho_1_minus_rho_2_plotting)),
                      #signal_name = paste0('rho_2 = ', as.vector(rho_2_plotting)),
                      ci_width = as.vector(ci_width),
                      gamma = as.vector(gamma_index_plotting))

legend_colors <- c('0.2' = 'gold3', '0.25' = 'darkseagreen4',
                   '0.3' = 'darkslategray4', '0.35' = 'deeppink4',
                   '0.4' = 'firebrick3')
# legend_colors = c('rho_2 = 20' = 'springgreen4', 'rho_2 = 25' = 'violetred4', 'rho_2 = 30' = 'slateblue4')
figure2a <- ggplot()
for (ri in 1:length(rho_2_check)) {

  # signal_name = paste0('rho2 = ', rho_2_check[ri])

  df_subset <- plot_df[plot_df$rho_2 == rho_2_check[ri], ]
  figure2a <- figure2a +
    geom_line(aes(x = gamma, y = ci_width, color = rho_1_minus_rho_2), size = 1.0, alpha = 0.7, data = df_subset)
}
figure2a <- figure2a +
  xlab(TeX('$\\gamma$')) + ylab('Average CI width') +
  scale_color_manual(breaks = c('0.2', '0.25', '0.3', '0.35', '0.4'), values = legend_colors) +
  # ggtitle(TeX('Average confidence interval width as a function of $\\gamma$')) +
  labs(color = TeX('$\\rho_1 - \\rho_2$')) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
figure2a
ggsave('figures/conf_width_gamma_bernoulli.pdf',
       plot = figure2a, device = 'pdf', width = 4.5, height = 3.5)

# Figure 2b - Tradeoff of signal regime, RAND index for clustering, separated by gamma
# --------------------------------------------------------------------------------------
#ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(xi_variances_thinning[n_index, K_true_index, K_index, , , ])
# ci_width <- apply(ci_width, c(1, 2), mean)
rand_slice <- apply(rand_results[n_index, K_true_index, K_index, , , ], c(1, 2), mean)
gamma_index_plotting <- matrix(rep(gamma_check, length(signal_regimes)), ncol = length(gamma_check), byrow = TRUE)
rho_2_plotting <- matrix(rep(rho_2_check, length(gamma_check)), ncol = length(gamma_check))
rho_1_minus_rho_2_plotting <- matrix(rep(rho_1_minus_rho_2, length(gamma_check)), ncol = length(gamma_check))
signal_index_plotting <- matrix(rep(1:length(signal_regimes), length(gamma_check)), ncol = length(gamma_check))

legend_colors <- c('0.2' = 'gold3', '0.25' = 'darkseagreen4',
                   '0.3' = 'darkslategray4', '0.35' = 'deeppink4',
                   '0.4' = 'firebrick3')

plot_df <- data.frame(signal_index = as.vector(signal_index_plotting),
                      rho_2 = as.character(as.vector(rho_2_plotting)),
                      rho_1_minus_rho_2 = as.character(as.vector(rho_1_minus_rho_2_plotting)),
                      rand_index = as.vector(rand_slice),
                      gamma = as.vector(gamma_index_plotting))
figure2b <- ggplot(plot_df) +
  geom_line(aes(x = gamma, y = rand_index, color = rho_1_minus_rho_2), size = 0.9) +
  scale_color_manual(breaks = c('100', '101', '102', '103', '104'), values = legend_colors) + # breaks are wrong because I don't want a legend to appear...
  xlab(TeX('$\\gamma')) + ylab('Adjusted RAND Index') +
  coord_fixed() +
  ylim(c(0, 1.0)) +
  #ggtitle(TeX('Classification accuracy (Adjusted RAND index) as a function of $\\gamma$')) +
  labs(color = TeX('$\\rho_1 - \\rho_2$')) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
figure2b
ggsave('figures/randindex_gamma_bernoulli.pdf',
       plot = figure2b, device = 'pdf', width = 4.5, height = 3.5)
