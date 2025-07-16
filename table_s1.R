# ================================================= #
# == R code to replicate the results in Table S1 == #
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
eps_check <- c(0.5)
signal_regimes <- list(c(30, 27)) # Signal regimes for mean matrix
tau2 <- 25
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

estimates_thinning_gaussian <- array(0, dim = c(length(n_check),
                                                length(K_true_check),
                                                length(K_check),
                                                length(signal_regimes),
                                                length(eps_check),
                                                num_sim))
targets_thinning_gaussian <- array(0, dim = c(length(n_check),
                                              length(K_true_check),
                                              length(K_check),
                                              length(signal_regimes),
                                              length(eps_check),
                                              num_sim))
variances_thinning_gaussian <- array(0, dim = c(length(n_check),
                                                length(K_true_check),
                                                length(K_check),
                                                length(signal_regimes),
                                                length(eps_check),
                                                num_sim))

estimates_naive_gaussian <- array(0, dim = c(length(n_check),
                                             length(K_true_check),
                                             length(K_check),
                                             length(signal_regimes),
                                             length(eps_check),
                                             num_sim))
targets_naive_gaussian <- array(0, dim = c(length(n_check),
                                           length(K_true_check),
                                           length(K_check),
                                           length(signal_regimes),
                                           length(eps_check),
                                           num_sim))
variances_naive_gaussian <- array(0, dim = c(length(n_check),
                                             length(K_true_check),
                                             length(K_check),
                                             length(signal_regimes),
                                             length(eps_check),
                                             num_sim))
covered_thinning_gaussian <- array(0, dim = c(length(n_check),
                                              length(K_true_check),
                                              length(K_check),
                                              length(signal_regimes),
                                              length(eps_check)))
covered_naive_gaussian <- array(0, dim = c(length(n_check),
                                           length(K_true_check),
                                           length(K_check),
                                           length(signal_regimes),
                                           length(eps_check)))

estimates_thinning_poisson <- array(0, dim = c(length(n_check),
                                               length(K_true_check),
                                               length(K_check),
                                               length(signal_regimes),
                                               length(eps_check),
                                               num_sim))
targets_thinning_poisson <- array(0, dim = c(length(n_check),
                                             length(K_true_check),
                                             length(K_check),
                                             length(signal_regimes),
                                             length(eps_check),
                                             num_sim))
variances_thinning_poisson <- array(0, dim = c(length(n_check),
                                               length(K_true_check),
                                               length(K_check),
                                               length(signal_regimes),
                                               length(eps_check),
                                               num_sim))

estimates_naive_poisson <- array(0, dim = c(length(n_check),
                                            length(K_true_check),
                                            length(K_check),
                                            length(signal_regimes),
                                            length(eps_check),
                                            num_sim))
targets_naive_poisson <- array(0, dim = c(length(n_check),
                                          length(K_true_check),
                                          length(K_check),
                                          length(signal_regimes),
                                          length(eps_check),
                                          num_sim))
variances_naive_poisson <- array(0, dim = c(length(n_check),
                                            length(K_true_check),
                                            length(K_check),
                                            length(signal_regimes),
                                            length(eps_check),
                                            num_sim))
covered_thinning_poisson <- array(0, dim = c(length(n_check),
                                             length(K_true_check),
                                             length(K_check),
                                             length(signal_regimes),
                                             length(eps_check)))
covered_naive_poisson <- array(0, dim = c(length(n_check),
                                          length(K_true_check),
                                          length(K_check),
                                          length(signal_regimes),
                                          length(eps_check)))

# ==============
# == Simulate ==
# ==============

cat("Simulating for Table S1\n")

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

        # Different values of epsilon in thinning
        for (eps_index in 1:length(eps_check)) {
          eps <- eps_check[eps_index]
          cat(paste0('\t\t\t\t(eps = ', eps, ')\n'))

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
            A_gaussian <- matrix(rnorm(n = n^2, mean = as.vector(M), sd = sqrt(tau2)), nrow = n)
            A_poisson <- matrix(rpois(n = n^2, lambda = as.vector(M)), nrow = n)

            # Split A
            A_split_gaussian <- networkinference::split_matrix(A_gaussian, distribution = 'gaussian',
                                                               epsilon = eps, allow_self_loops = TRUE,
                                                               is_directed = TRUE, tau = sqrt(tau2))
            A_split_poisson <- networkinference::split_matrix(A_poisson, distribution = "poisson",
                                                      epsilon = eps)


            Atr_gaussian <- A_split_gaussian$Atr
            Ate_gaussian <- A_split_gaussian$Ate
            Atr_poisson <- A_split_poisson$Atr
            Ate_poisson <- A_split_poisson$Ate

            # Clustering for naive and thinning
            z_hat_gaussian <- nett::spec_clust(Atr_gaussian, K = K)
            z_hat_naive_gaussian <- nett::spec_clust(A_gaussian, K = K)
            z_hat_poisson <- nett::spec_clust(Atr_poisson, K = K)
            z_hat_naive_poisson <- nett::spec_clust(A_poisson, K = K)

            # Estimator, estimator variance, and target for the proposed approach
            thinning_infer_gaussian <-
              networkinference::infer_network(Ate = Ate_gaussian, u = u,
                                              communities = z_hat_gaussian,
                                              distribution = "gaussian",
                                              K = K, epsilon = eps, tau = sqrt(tau2))
            thinning_infer_poisson <-
              networkinference::infer_network(Ate = Ate_poisson, u = u,
                                              communities = z_hat_poisson,
                                              distribution = "poisson",
                                              K = K, epsilon = eps, tau = sqrt(tau2))

            estimate_gaussian <- thinning_infer_gaussian$estimate
            estimate_poisson <- thinning_infer_poisson$estimate

            estimate_var_gaussian <- thinning_infer_gaussian$estimate_variance
            estimate_var_poisson <- thinning_infer_poisson$estimate_variance

            target_gaussian <- networkinference::check_target_of_inference(M = M, u = u,
                                                                           communities = z_hat_gaussian,
                                                                           K = K)
            target_poisson <- networkinference::check_target_of_inference(M = M, u = u,
                                                                          communities = z_hat_poisson,
                                                                          K = K)

            # Estimator, estimator variance, and target for the naive approach
            naive_infer_gaussian <-
              networkinference::infer_network(Ate = A_gaussian, u = u,
                                              communities = z_hat_naive_gaussian,
                                              distribution = "gaussian",
                                              K = K, epsilon = 0, tau = sqrt(tau2))
            naive_infer_poisson <-
              networkinference::infer_network(Ate = A_poisson, u = u,
                                              communities = z_hat_naive_poisson,
                                              distribution = "poisson",
                                              K = K, epsilon = 0)

            estimate_naive_gaussian <- naive_infer_gaussian$estimate
            estimate_naive_poisson <- naive_infer_poisson$estimate

            estimate_var_naive_gaussian <- naive_infer_gaussian$estimate_variance
            estimate_var_naive_poisson <- naive_infer_poisson$estimate_variance

            target_naive_gaussian <- networkinference::check_target_of_inference(M = M, u = u,
                                                                        communities = z_hat_naive_gaussian,
                                                                        K = K)
            target_naive_poisson <- networkinference::check_target_of_inference(M = M, u = u,
                                                                                 communities = z_hat_naive_poisson,
                                                                                 K = K)

            # Save results
            # ------------
            estimates_thinning_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_gaussian
            targets_thinning_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target_gaussian
            variances_thinning_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var_gaussian
            estimates_naive_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_naive_gaussian
            targets_naive_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target_naive_gaussian
            variances_naive_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var_naive_gaussian

            estimates_thinning_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_poisson
            targets_thinning_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target_poisson
            variances_thinning_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var_poisson
            estimates_naive_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_naive_poisson
            targets_naive_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target_naive_poisson
            variances_naive_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var_naive_poisson

            # Do inference and save results
            margin_of_error_gaussian <- qnorm(1 - alpha / 2) * sqrt(estimate_var_gaussian)
            margin_of_error_naive_gaussian <- qnorm(1 - alpha / 2) * sqrt(estimate_var_naive_gaussian)
            margin_of_error_poisson <- qnorm(1 - alpha / 2) * sqrt(estimate_var_poisson)
            margin_of_error_naive_poisson <- qnorm(1 - alpha / 2) * sqrt(estimate_var_naive_poisson)

            # Check coverage (thinning and naive)
            covered_thinning_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_thinning_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index] +
              as.integer((estimate_gaussian - margin_of_error_gaussian <= target_gaussian) & (target_gaussian <= estimate_gaussian + margin_of_error_gaussian))
            covered_thinning_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_thinning_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index] +
              as.integer((estimate_poisson - margin_of_error_poisson <= target_poisson) & (target_poisson <= estimate_poisson + margin_of_error_poisson))

            covered_naive_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_naive_gaussian[n_index, K_true_index, K_index, signal_regime_index, eps_index] +
              as.integer((estimate_naive_gaussian - margin_of_error_naive_gaussian <= target_naive_gaussian) & (target_naive_gaussian <= estimate_naive_gaussian + margin_of_error_naive_gaussian))
            covered_naive_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_naive_poisson[n_index, K_true_index, K_index, signal_regime_index, eps_index] +
              as.integer((estimate_naive_poisson - margin_of_error_naive_poisson <= target_naive_poisson) & (target_naive_poisson <= estimate_naive_poisson + margin_of_error_naive_poisson))
          }
        }
      }
    }
  }
}
covered_thinning_gaussian <- covered_thinning_gaussian / num_sim
covered_naive_gaussian <- covered_naive_gaussian / num_sim
covered_thinning_poisson <- covered_thinning_poisson / num_sim
covered_naive_poisson <- covered_naive_poisson / num_sim

# ================================================ #
# == Save/load results (saves computation time) == #
# ================================================ #

# saveRDS(estimates_thinning_gaussian, file = "saved_simulation_data/table_s1_estimates_thinning_gaussian.RDS")
# saveRDS(targets_thinning_gaussian, file = "saved_simulation_data/table_s1_targets_thinning_gaussian.RDS")
# saveRDS(variances_thinning_gaussian, file = "saved_simulation_data/table_s1_variances_thinning_gaussian.RDS")
# saveRDS(estimates_naive_gaussian, file = "saved_simulation_data/table_s1_estimates_naive_gaussian.RDS")
# saveRDS(targets_naive_gaussian, file = "saved_simulation_data/table_s1_targets_naive_gaussian.RDS")
# saveRDS(variances_naive_gaussian, file = "saved_simulation_data/table_s1_variances_naive_gaussian.RDS")
# saveRDS(estimates_thinning_poisson, file = "saved_simulation_data/table_s1_estimates_thinning_poisson.RDS")
# saveRDS(targets_thinning_poisson, file = "saved_simulation_data/table_s1_targets_thinning_poisson.RDS")
# saveRDS(variances_thinning_poisson, file = "saved_simulation_data/table_s1_variances_thinning_poisson.RDS")
# saveRDS(estimates_naive_poisson, file = "saved_simulation_data/table_s1_estimates_naive_poisson.RDS")
# saveRDS(targets_naive_poisson, file = "saved_simulation_data/table_s1_targets_naive_poisson.RDS")
# saveRDS(variances_naive_poisson, file = "saved_simulation_data/table_s1_variances_naive_poisson.RDS")
# saveRDS(covered_thinning_gaussian, file = "saved_simulation_data/table_s1_covered_thinning_gaussian.RDS")
# saveRDS(covered_naive_gaussian, file = "saved_simulation_data/table_s1_covered_naive_gaussian.RDS")
# saveRDS(covered_thinning_poisson, file = "saved_simulation_data/table_s1_covered_thinning_poisson.RDS")
# saveRDS(covered_naive_poisson, file = "saved_simulation_data/table_s1_covered_naive_poisson.RDS")

estimates_thinning_gaussian <- readRDS("saved_simulation_data/table_s1_estimates_thinning_gaussian.RDS")
targets_thinning_gaussian <- readRDS("saved_simulation_data/table_s1_targets_thinning_gaussian.RDS")
variances_thinning_gaussian <- readRDS("saved_simulation_data/table_s1_variances_thinning_gaussian.RDS")
estimates_naive_gaussian <- readRDS("saved_simulation_data/table_s1_estimates_naive_gaussian.RDS")
targets_naive_gaussian <- readRDS("saved_simulation_data/table_s1_targets_naive_gaussian.RDS")
variances_naive_gaussian <- readRDS("saved_simulation_data/table_s1_variances_naive_gaussian.RDS")
estimates_thinning_poisson <- readRDS("saved_simulation_data/table_s1_estimates_thinning_poisson.RDS")
targets_thinning_poisson <- readRDS("saved_simulation_data/table_s1_targets_thinning_poisson.RDS")
variances_thinning_poisson <- readRDS("saved_simulation_data/table_s1_variances_thinning_poisson.RDS")
estimates_naive_poisson <- readRDS("saved_simulation_data/table_s1_estimates_naive_poisson.RDS")
targets_naive_poisson <- readRDS("saved_simulation_data/table_s1_targets_naive_poisson.RDS")
variances_naive_poisson <- readRDS("saved_simulation_data/table_s1_variances_naive_poisson.RDS")
covered_thinning_gaussian <- readRDS("saved_simulation_data/table_s1_covered_thinning_gaussian.RDS")
covered_naive_gaussian <- readRDS("saved_simulation_data/table_s1_covered_naive_gaussian.RDS")
covered_thinning_poisson <- readRDS("saved_simulation_data/table_s1_covered_thinning_poisson.RDS")
covered_naive_poisson <- readRDS("saved_simulation_data/table_s1_covered_naive_poisson.RDS")

# ================= #
# == Print table == #
# ================= #

plot_table_s1 <- function(data) {

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

cat('GAUSSIAN\n')
cat('--------\n')
cat('Naive\n')
plot_table_s1(covered_naive_gaussian)

cat('Proposed\n')
plot_table_s1(covered_thinning_gaussian)

cat('\n')
cat('POISSON\n')
cat('--------\n')
cat('Naive\n')
plot_table_s1(covered_naive_poisson)

cat('Proposed\n')
plot_table_s1(covered_thinning_poisson)
