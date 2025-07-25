# ================================================= #
# == R code to replicate the results in Table S2 == #
# == Author: Ethan Ancell                        == #
# ================================================= #

library("networkinference")
library("tidyverse")
library("nett")
library("mclust")
library("latex2exp")
library("networkinference")

ggplot2::theme_set(theme_minimal())

# =========================
# == Simulation settings ==
# =========================

set.seed(1)
num_sim <- 5000 # Number of replications
n_check <- c(100, 200, 500) # Total network size
K_check <- c(2, 5, 10) # Guessed number of communities
K_true_check <- c(5)
gamma_check <- c(0.25)
signal_regimes <- list(c(0.75, 0.5)) # Signal regimes for mean matrix
alpha <- 0.10

# ====================== #
# == Simulation setup == #
# ====================== #

# Unpack rho_1 and rho_2 from the signal regimes
rho_1 <- signal_regimes[[1]][1] # Assumes that this stays the same across all of them
rho_2_check <- rep(NA, length(signal_regimes))
for (i in 1:length(signal_regimes)) {
  rho_2_check[i] <- signal_regimes[[i]][2]
}
rho_1_minus_rho_2 <- rho_1 - rho_2_check

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

cat("Simulating for Table S2\n")

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

            # Linear combination vector
            u <- c(1, rep(0, K^2 - 1))

            # Draw A
            A <- matrix(rbinom(n = n^2, size = 1, prob = as.vector(M)), nrow = n)

            # Split A
            A_split <- networkinference::split_network(A, distribution = "bernoulli",
                                                      gamma = gamma)

            Atr <- A_split$Atr
            Ate <- A_split$Ate

            # Clustering for naive and thinning
            z_hat <- nett::spec_clust(Atr, K = K)
            z_hat_naive <- nett::spec_clust(A, K = K)

            # Estimator, estimator variance, and target for the proposed approach
            fission_infer <-
              networkinference::infer_network(Ate = Ate, u = u,
                                              communities = z_hat,
                                              distribution = "bernoulli",
                                              gamma = gamma,
                                              Atr = Atr)

            estimate <- fission_infer$estimate
            estimate_var <- fission_infer$estimate_variance
            target_theta <- networkinference::check_target_of_inference(M = M, u = u,
                                                                        communities = z_hat)
            target_xi <- networkinference::check_target_of_inference(M = M, u = u,
                                                                     communities = z_hat,
                                                                     bernoulli_target = TRUE,
                                                                     gamma = gamma,
                                                                     Atr = Atr)

            # Estimator, estimator variance, and target for the naive approach
            # (Note that this is manually implemented, as my R package
            # does not implement this functionality.)

            # Helper matrices
            z_hat_naive_matrix <- matrix(rep(NA, n*K), nrow = n)
            for (i in 1:K) {
              z_hat_naive_matrix[, i] <- 1 * (z_hat_naive == i)
            }
            n_hat_naive <- apply(z_hat_naive_matrix, 2, sum)
            NN_inv_naive <- diag(1 / diag(t(z_hat_naive_matrix) %*% z_hat_naive_matrix))
            comm_pair_sample_size_naive <- t(z_hat_naive_matrix) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_naive_matrix

            # Estimate and estimator
            theta_naive_matrix <- NN_inv_naive %*% t(z_hat_naive_matrix) %*% M %*% z_hat_naive_matrix %*% NN_inv_naive
            theta_naive_ <- t(u) %*% as.vector(theta_naive_matrix)
            theta_est_naive_matrix <- NN_inv_naive %*% t(z_hat_naive_matrix) %*% A %*% z_hat_naive_matrix %*% NN_inv_naive
            theta_est_naive_ <- t(u) %*% as.vector(theta_est_naive_matrix)

            # Estimator variance
            theta_var_naive_matrix <- (theta_est_naive_matrix * (1 - theta_est_naive_matrix)) / comm_pair_sample_size_naive
            theta_var_naive_ <- t(u) %*% diag(as.vector(theta_var_naive_matrix)) %*% u

            # Do inference and save results
            margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)
            margin_of_error_naive <- qnorm(1 - alpha / 2) * sqrt(theta_var_naive_)

            # Check and save coverage
            xi_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] <-
              xi_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] +
              as.integer((estimate - margin_of_error <= target_xi) & (target_xi <= estimate + margin_of_error))
            theta_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] <-
              theta_covered_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index] +
              as.integer((estimate - margin_of_error <= target_theta) & (target_theta <= estimate + margin_of_error))
            theta_covered_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index] <-
              theta_covered_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index] +
              as.integer((theta_est_naive_ - margin_of_error_naive <= theta_naive_) & (theta_naive_ <= theta_est_naive_ + margin_of_error_naive))
          }
        }
      }
    }
  }
}
xi_covered_thinning <- xi_covered_thinning / num_sim
theta_covered_thinning <- theta_covered_thinning / num_sim
theta_covered_naive <- theta_covered_naive / num_sim

# ================================================ #
# == Save/load results (saves computation time) == #
# ================================================ #

# saveRDS(xi_covered_thinning, file = "saved_simulation_data/table_s1_xi_covered_thinning.RDS")
# saveRDS(theta_covered_thinning, file = "saved_simulation_data/table_s1_theta_covered_thinning.RDS")
# saveRDS(theta_covered_naive, file = "saved_simulation_data/table_s1_theta_covered_naive.RDS")

xi_covered_thinning <- readRDS("saved_simulation_data/table_s1_xi_covered_thinning.RDS")
theta_covered_thinning <- readRDS("saved_simulation_data/table_s1_theta_covered_thinning.RDS")
theta_covered_naive <- readRDS("saved_simulation_data/table_s1_theta_covered_naive.RDS")

# ================= #
# == Print table == #
# ================= #

plot_table_s2 <- function(data) {

  signal_index <- 1
  K_true_index <- 1
  eps_index <- 1

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

      cat(sprintf("   %.2f |", 100 * data[n_index, K_true_index, K_index, signal_index, eps_index]))
    }
    cat('\n')
  }
}

cat('Bernoulli - Naive\n')
plot_table_s2(theta_covered_naive)

cat('Bernoulli - Xi - Thinning\n')
plot_table_s2(xi_covered_thinning)

cat('Bernoulli - Theta - Thinning\n')
plot_table_s2(theta_covered_thinning)
