# ============================================
# == Simulations for network thinning paper ==
# == (Poisson edges)                        ==
# == Author: Ethan Ancell                   ==
# ============================================

# Libraries
library(ggplot2)
library(tidyverse)
library(nett)
library(mclust)
library(latex2exp)

ggplot2::theme_set(theme_minimal())

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
# eps_check <- c(0.2, 0.5, 0.7)


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
eps_check <- c(0.5)
signal_regimes <- list(c(30, 27)) # Signal regimes for mean matrix

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
eps_check <- c(0.5)
signal_regimes <- list(c(30, 27)) # Signal regimes for mean matrix

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
eps_check <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
signal_regimes <- list(c(30, 29), c(30, 27), c(30, 25), c(30, 23), c(30, 21)) # Signal regimes for mean matrix

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
estimates_thinning <- array(0, dim = c(length(n_check),
                                       length(K_true_check),
                                       length(K_check),
                                       length(signal_regimes),
                                       length(eps_check),
                                       num_sim))
targets_thinning <- array(0, dim = c(length(n_check),
                                     length(K_true_check),
                                     length(K_check),
                                     length(signal_regimes),
                                     length(eps_check),
                                     num_sim))
variances_thinning <- array(0, dim = c(length(n_check),
                                       length(K_true_check),
                                       length(K_check),
                                       length(signal_regimes),
                                       length(eps_check),
                                       num_sim))

estimates_naive <- array(0, dim = c(length(n_check),
                                    length(K_true_check),
                                    length(K_check),
                                    length(signal_regimes),
                                    length(eps_check),
                                    num_sim))
targets_naive <- array(0, dim = c(length(n_check),
                                  length(K_true_check),
                                  length(K_check),
                                  length(signal_regimes),
                                  length(eps_check),
                                  num_sim))
variances_naive <- array(0, dim = c(length(n_check),
                                    length(K_true_check),
                                    length(K_check),
                                    length(signal_regimes),
                                    length(eps_check),
                                    num_sim))
rand_results <- array(0, dim = c(length(n_check),
                                 length(K_true_check),
                                 length(K_check),
                                 length(signal_regimes),
                                 length(eps_check),
                                 num_sim))

# Coverage results of (1-alpha) confidence intervals
covered_thinning <- array(0, dim = c(length(n_check),
                                     length(K_true_check),
                                     length(K_check),
                                     length(signal_regimes),
                                     length(eps_check)))
covered_naive <- array(0, dim = c(length(n_check),
                                  length(K_true_check),
                                  length(K_check),
                                  length(signal_regimes),
                                  length(eps_check)))

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

            # Sample u
            if (use_random_u) {
              # Sample the vector u from the unit K^2-sphere
              u <- rnorm(K^2)
              u <- u / sqrt(sum(u^2))
            } else {
              u <- c(1, rep(0, K^2 - 1))
            }

            # Draw A
            A <- matrix(rpois(n = n^2, lambda = as.vector(M)), nrow = n)

            # Thin A
            A_tr <- matrix(rbinom(n = n^2, size = as.vector(A), prob = eps), nrow = n)
            A_te <- A - A_tr

            # Cluster (thinning)
            z_hat_initial <- nett::spec_clust(A_tr, K = K)
            z_hat <- matrix(rep(NA, n*K), nrow = n)
            for (i in 1:K) {
              z_hat[, i] <- 1 * (z_hat_initial == i)
            }
            n_hat <- apply(z_hat, 2, sum)

            # Check RAND index agreement between true clusters and non true clusters
            rand_results[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <- mclust::adjustedRandIndex(z_hat_initial, z_true)
            #rand_results[]

            # Cluster (naive)
            z_hat_initial_naive <- nett::spec_clust(A, K = K)
            z_hat_naive <- matrix(rep(NA, n*K), nrow = n)
            for (i in 1:K) {
              z_hat_naive[, i] <- 1 * (z_hat_initial_naive == i)
            }
            n_hat_naive <- apply(z_hat_naive, 2, sum)

            # Target and estimator (thinning)
            NN_inv <- diag(1 / diag(t(z_hat) %*% z_hat))
            estimate_blocky <- NN_inv %*% t(z_hat) %*% A_te %*% z_hat %*% NN_inv
            estimate <- (1 - eps)^(-1) * t(u) %*% as.vector(NN_inv %*% t(z_hat) %*% A_te %*% z_hat %*% NN_inv)
            target <- t(u) %*% as.vector(NN_inv %*% t(z_hat) %*% M %*% z_hat %*% NN_inv)

            # Target and estimator (naive)
            NN_inv_naive <- diag(1 / diag(t(z_hat_naive) %*% z_hat_naive))
            estimate_blocky_naive <- NN_inv_naive %*% t(z_hat_naive) %*% A %*% z_hat_naive %*% NN_inv_naive
            estimate_naive <- t(u) %*% as.vector(NN_inv_naive %*% t(z_hat_naive) %*% A %*% z_hat_naive %*% NN_inv_naive)
            target_naive <- t(u) %*% as.vector(NN_inv_naive %*% t(z_hat_naive) %*% M %*% z_hat_naive %*% NN_inv_naive)

            # Estimated variance of estimator (thinning and naive)
            Delta_n <- matrix(rep(NA, K^2), nrow = K)
            Delta_n_naive <- matrix(rep(NA, K^2), nrow = K)
            for (k in 1:K) {
              for (l in 1:K) {
                # Thinning
                k_nodes <- which(z_hat_initial == k)
                l_nodes <- which(z_hat_initial == l)

                Jn <- sum(A_te[k_nodes, l_nodes]) / estimate_blocky[k, l]^2
                # Vn <- sum((A_te[k_nodes, l_nodes] / estimate_blocky[k, l] - 1)^2)

                # Delta_n[k, l] <- Vn / Jn^2
                Delta_n[k, l] <- 1 / Jn

                # Naive
                k_nodes_naive <- which(z_hat_initial_naive == k)
                l_nodes_naive <- which(z_hat_initial_naive == l)

                Jn_naive <- sum(A[k_nodes_naive, l_nodes_naive]) / estimate_blocky_naive[k, l]^2
                # Vn_naive <- sum((A[k_nodes_naive, l_nodes_naive] / estimate_blocky_naive[k, l] - 1)^2)

                # Delta_n_naive[k, l] <- Vn_naive / Jn_naive^2
                Delta_n_naive[k, l] <- 1 / Jn_naive
              }
            }
            estimate_var <- (1-eps)^(-2) * t(u) %*% diag(as.vector(Delta_n)) %*% u
            estimate_var_naive <- t(u) %*% diag(as.vector(Delta_n_naive)) %*% u

            # Construct a confidence interval for the target of inference
            margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)
            margin_of_error_naive <- qnorm(1 - alpha / 2) * sqrt(estimate_var_naive)

            # Check coverage (thinning and naive)
            covered_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index] +
              as.integer((estimate - margin_of_error <= target) & (target <= estimate + margin_of_error))

            covered_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index] +
              as.integer((estimate_naive - margin_of_error_naive <= target_naive) & (target_naive <= estimate_naive + margin_of_error_naive))

            # Save everything else
            # --------------------
            estimates_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate
            targets_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target
            variances_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var
            estimates_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_naive
            targets_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target_naive
            variances_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var_naive
          }
        }
      }
    }
  }
}

covered_thinning <- covered_thinning / num_sim
covered_naive <- covered_naive / num_sim

covered_thinning
covered_naive

# ====================================================
# == SAVE the simulation results after running them ==
# == Uncomment/comment the lines appropriately.     ==
# ====================================================

# -----------
# - Table 1 -
# -----------
# saveRDS(estimates_thinning, file = 'save_data/SIM1_POISSON_estimates_thinning.RDS')
# saveRDS(targets_thinning, file = 'save_data/SIM1_POISSON_targets_thinning.RDS')
# saveRDS(variances_thinning, file = 'save_data/SIM1_POISSON_variances_thinning.RDS')
# saveRDS(estimates_naive, file = 'save_data/SIM1_POISSON_estimates_naive.RDS')
# saveRDS(targets_naive, file = 'save_data/SIM1_POISSON_targets_naive.RDS')
# saveRDS(variances_naive, file = 'save_data/SIM1_POISSON_variances_naive.RDS')
# saveRDS(rand_results, file = 'save_data/SIM1_POISSON_rand_results.RDS')
# saveRDS(covered_thinning, file = 'save_data/SIM1_POISSON_covered_thinning.RDS')
# saveRDS(covered_naive, file = 'save_data/SIM1_POISSON_covered_naive.RDS')

estimates_thinning <- readRDS('save_data/SIM1_POISSON_estimates_thinning.RDS')
targets_thinning <- readRDS('save_data/SIM1_POISSON_targets_thinning.RDS')
variances_thinning <- readRDS('save_data/SIM1_POISSON_variances_thinning.RDS')
estimates_naive <- readRDS('save_data/SIM1_POISSON_estimates_naive.RDS')
targets_naive <- readRDS('save_data/SIM1_POISSON_targets_naive.RDS')
variances_naive <- readRDS('save_data/SIM1_POISSON_variances_naive.RDS')
rand_results <- readRDS('save_data/SIM1_POISSON_rand_results.RDS')
covered_thinning <- readRDS('save_data/SIM1_POISSON_covered_thinning.RDS')
covered_naive <- readRDS('save_data/SIM1_POISSON_covered_naive.RDS')

# ------------
# - Figure 1 -
# ------------
# saveRDS(estimates_thinning, file = 'save_data/SIM2_POISSON_estimates_thinning.RDS')
# saveRDS(targets_thinning, file = 'save_data/SIM2_POISSON_targets_thinning.RDS')
# saveRDS(variances_thinning, file = 'save_data/SIM2_POISSON_variances_thinning.RDS')
# saveRDS(estimates_naive, file = 'save_data/SIM2_POISSON_estimates_naive.RDS')
# saveRDS(targets_naive, file = 'save_data/SIM2_POISSON_targets_naive.RDS')
# saveRDS(variances_naive, file = 'save_data/SIM2_POISSON_variances_naive.RDS')
# saveRDS(rand_results, file = 'save_data/SIM2_POISSON_rand_results.RDS')
# saveRDS(covered_thinning, file = 'save_data/SIM2_POISSON_covered_thinning.RDS')
# saveRDS(covered_naive, file = 'save_data/SIM2_POISSON_covered_naive.RDS')

estimates_thinning <- readRDS('save_data/SIM2_POISSON_estimates_thinning.RDS')
targets_thinning <- readRDS('save_data/SIM2_POISSON_targets_thinning.RDS')
variances_thinning <- readRDS('save_data/SIM2_POISSON_variances_thinning.RDS')
estimates_naive <- readRDS('save_data/SIM2_POISSON_estimates_naive.RDS')
targets_naive <- readRDS('save_data/SIM2_POISSON_targets_naive.RDS')
variances_naive <- readRDS('save_data/SIM2_POISSON_variances_naive.RDS')
rand_results <- readRDS('save_data/SIM2_POISSON_rand_results.RDS')
covered_thinning <- readRDS('save_data/SIM2_POISSON_covered_thinning.RDS')
covered_naive <- readRDS('save_data/SIM2_POISSON_covered_naive.RDS')

# ------------
# - Figure 2 -
# ------------
# saveRDS(estimates_thinning, file = 'save_data/SIM3_POISSON_estimates_thinning.RDS')
# saveRDS(targets_thinning, file = 'save_data/SIM3_POISSON_targets_thinning.RDS')
# saveRDS(variances_thinning, file = 'save_data/SIM3_POISSON_variances_thinning.RDS')
# saveRDS(estimates_naive, file = 'save_data/SIM3_POISSON_estimates_naive.RDS')
# saveRDS(targets_naive, file = 'save_data/SIM3_POISSON_targets_naive.RDS')
# saveRDS(variances_naive, file = 'save_data/SIM3_POISSON_variances_naive.RDS')
# saveRDS(rand_results, file = 'save_data/SIM3_POISSON_rand_results.RDS')
# saveRDS(covered_thinning, file = 'save_data/SIM3_POISSON_covered_thinning.RDS')
# saveRDS(covered_naive, file = 'save_data/SIM3_POISSON_covered_naive.RDS')

estimates_thinning <- readRDS('save_data/SIM3_POISSON_estimates_thinning.RDS')
targets_thinning <- readRDS('save_data/SIM3_POISSON_targets_thinning.RDS')
variances_thinning <- readRDS('save_data/SIM3_POISSON_variances_thinning.RDS')
estimates_naive <- readRDS('save_data/SIM3_POISSON_estimates_naive.RDS')
targets_naive <- readRDS('save_data/SIM3_POISSON_targets_naive.RDS')
variances_naive <- readRDS('save_data/SIM3_POISSON_variances_naive.RDS')
rand_results <- readRDS('save_data/SIM3_POISSON_rand_results.RDS')
covered_thinning <- readRDS('save_data/SIM3_POISSON_covered_thinning.RDS')
covered_naive <- readRDS('save_data/SIM3_POISSON_covered_naive.RDS')

# ==================================
# == Plotting and viewing results ==
# ==================================

# Table 1 - Checking coverage of thinning and naive
# -------------------------------------------------
plot_table_1 <- function(data) {

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

cat('Naive\n')
plot_table_1(covered_naive)

cat('Thinning\n')
plot_table_1(covered_thinning)

# Figure 1 - Showing a comparison of nominal vs empirical coverage
# ----------------------------------------------------------------

checked_coverages <- seq(0, 1.0, length.out = 100)
n_coverages <- length(checked_coverages)

n_index <- 1
K_true_index <- 1
K_true <- K_true_check[K_true_index]
signal_index <- 1

# Results
figure1_coverage_thinning <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))
figure1_coverage_naive <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))

for (K_index in 1:length(K_check)) {

  for (eps_index in 1:length(eps_check)) {

    for (coverage_index in 1:length(checked_coverages)) {

      coverage <- checked_coverages[coverage_index]
      z_quantile <- qnorm(1 - (1 - coverage) / 2)

      # Extract relevant slices from the results array
      targets_thinning_slice <- targets_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ]
      targets_naive_slice <- targets_naive[n_index, K_true_index, K_index, signal_index, eps_index, ]

      upper_bounds_thinning <-
        estimates_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ] +
        z_quantile * sqrt(variances_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ])
      lower_bounds_thinning <-
        estimates_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ] -
        z_quantile * sqrt(variances_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ])

      upper_bounds_naive <-
        estimates_naive[n_index, K_true_index, K_index, signal_index, eps_index, ] +
        z_quantile * sqrt(variances_naive[n_index, K_true_index, K_index, signal_index, eps_index, ])
      lower_bounds_naive <-
        estimates_naive[n_index, K_true_index, K_index, signal_index, eps_index, ] -
        z_quantile * sqrt(variances_naive[n_index, K_true_index, K_index, signal_index, eps_index, ])

      # Store coverage
      figure1_coverage_thinning[K_index, eps_index, coverage_index] <-
        mean((lower_bounds_thinning <= targets_thinning_slice) & (targets_thinning_slice <= upper_bounds_thinning))
      figure1_coverage_naive[K_index, eps_index, coverage_index] <-
        mean((lower_bounds_naive <= targets_naive_slice) & (targets_naive_slice <= upper_bounds_naive))

    }

  }

}

# Stack all this information in a data frame that can be used for plotting
plot_df <- data.frame()
for (K_index in 1:length(K_check)) {
  for (eps_index in 1:length(eps_check)) {

    # Create a method name for thinning / naive that includes K = ---
    # thinning_method_name <- paste0('Thinning (K = ', K_check[K_index], ')')
    # naive_method_name <- paste0('Naive (K = ', K_check[K_index], ')')

    # Thinning
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_thinning[K_index, eps_index, ],
                                method = rep('Proposed', n_coverages),
                                # method_print = rep(thinning_method_name, n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages)))

    # Naive
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_naive[K_index, eps_index, ],
                                method = rep('Naive', n_coverages),
                                # method_print = rep(naive_method_name, n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages)))
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
ggsave('figures/nominal_empirical_coverage_poisson.pdf',
       plot = figure1, device = 'pdf', width = 4.5, height = 3.5)
#
# # Create plot
# # plot_df_thinning <- plot_df %>% filter(method == 'Proposed')
# # plot_df_naive <- plot_df %>% filter(method == 'Naive')
# legend_colors <- c('2' = 'steelblue1', '5' = 'steelblue3', '10' = 'steelblue4')
#
# figure1 <- ggplot(plot_df_thinning) +
#   geom_line(aes(x = nominal_coverage, y = coverage, color = K), linetype = 'solid', size = 0.6)
# figure1
#
# # Create plot
# # plot_df <- data.frame(index = 1:length(checked_coverages),
# #                       nominal_coverage = checked_coverages,
# #                       Naive = figure1_coverage_naive,
# #                       Proposed = figure1_coverage_thinning,
# #                       Ideal = checked_coverages)
# # plot_df <- plot_df %>% pivot_longer(c('Naive', 'Proposed', 'Ideal'))
#
# # plot_df_thinning <- plot_df %>% filter(method == 'Proposed')
# # plot_df_naive <- plot_df %>% filter(method == 'Naive', eps == 0.5)
# # ideal_line <- data.frame(nominal_coverage = checked_coverages, coverage = checked_coverages)
#
# # theme_set(theme_minimal())
# # legend_colors <- c('Ideal' = 'black', 'Naive' = 'skyblue', Proposed = 'coral')
# # legend_colors <- c('Ideal' = 'black',
# #                    'Naive (K = 2)' = 'sienna1', 'Naive (K = 5)' = 'sienna3', 'Naive (K = 10)' = 'sienna4',
# #                    'Proposed (K = 2)' = 'steelblue1', 'Proposed (K = 5)' = 'steelblue3', 'Proposed (K = 10)' = 'steelblue4')
# legend_colors <- c('2' = 'steelblue1', '5' = 'steelblue3', '10' = 'steelblue4')
#
# figure1 <- ggplot()
# for (K_index in 1:length(K_check)) {
#   K <- K_check[K_index]
#
#   # for (eps in eps_check) {
#   #   # Subset data
#   #   df_subset <- plot_df_thinning[plot_df_thinning$eps == eps & plot_df_thinning$K == K, ]
#   #
#   #   # Add thinning
#   #   figure1 <- figure1 +
#   #     geom_line(aes(x = nominal_coverage, y = coverage, color = method_print, linetype = as.character(eps)), size = 0.6, alpha = 0.8, data = df_subset)
#   # }
#
#   df_subset <- plot_df[plot_df]
#
#   # Subset data
#   df_subset <- plot_df_thinning[plot_df_thinning$K == K, ]
#
#   # Add thinning
#   figure1 <- figure1 +
#     geom_line(aes(x = nominal_coverage, y = coverage, color = method_print, linetype = as.character(eps)), size = 0.6, alpha = 0.8, data = df_subset)
#
#   # Add naive
#   df_subset <- plot_df_naive[plot_df_naive$K == K, ]
#
#   # Add thinning
#   figure1 <- figure1 +
#     geom_line(aes(x = nominal_coverage, y = coverage, color = method_print), size = 0.6, data = df_subset)
# }
#
# figure1 <- figure1 +
#   geom_line(aes(x = nominal_coverage, y = coverage, color = 'Ideal'), size = 0.4, data = ideal_line) +
#   scale_color_manual(breaks = c('Ideal', 'Naive (K = 2)', 'Naive (K = 5)',
#                                 'Naive (K = 10)', 'Proposed (K = 2)',
#                                 'Proposed (K = 5)', 'Proposed (K = 10)'),
#                      values = legend_colors) +
#   xlab('Nominal Coverage') + ylab('Empirical Coverage') +
#   coord_fixed() +
#   labs(color = 'Method', linetype = "Epsilon") +
#   theme(plot.title = element_text(hjust = 0.5))
#
# figure1


# ggplot(plot_df) +
#   geom_line(aes(x = nominal_coverage, y = coverage, color = 'Proposed'), data = plot_df_thinning) +
#   geom_line(aes(x = nominal_coverage, y = coverage, color = 'Naive'), data = plot_df_naive) +
#   xlab('Nominal Coverage') + ylab('Empirical Coverage') +
#   ggtitle(paste0('n = ', n_check[n_index], ', K_true = ', K_true, ', K = ',
#                  K_check[K_index], ', eps = ', eps_check[eps_index])) +
#   labs(color = 'Method')



# Figure 2 - Tradeoff of signal regime, epsilon, RAND index, CI width
# -------------------------------------------------------------------

n_index <- 1
K_true_index <- 1
K_index <- 1
alpha <- 0.10

# Figure 2a - CI width as a function of signal, separate by epsilon
# -----------------------------------------------------------------
ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_true_index, K_index, , , ])
# ci_width <- (2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_true_index, K_index, , , ])) / sqrt(estimates_thinning[n_index, K_true_index, K_index, , , ])
ci_width <- apply(ci_width, c(1, 2), mean)
eps_index_plotting <- matrix(rep(eps_check, length(signal_regimes)), ncol = length(eps_check), byrow = TRUE)
rho_2_plotting <- matrix(rep(rho_2_check, length(eps_check)), ncol = length(eps_check))
rho_1_minus_rho_2_plotting <- matrix(rep(rho_1_minus_rho_2, length(eps_check)), ncol = length(eps_check))
signal_index_plotting <- matrix(rep(1:length(signal_regimes), length(eps_check)), ncol = length(eps_check))

plot_df <- data.frame(signal_index = as.vector(signal_index_plotting),
                      rho_2 = as.vector(rho_2_plotting),
                      rho_1_minus_rho_2 = as.character(as.vector(rho_1_minus_rho_2_plotting)),
                      #signal_name = paste0('rho_2 = ', as.vector(rho_2_plotting)),
                      ci_width = as.vector(ci_width),
                      eps = as.vector(eps_index_plotting))

legend_colors <- c('1' = 'gold3', '3' = 'darkseagreen4',
                   '5' = 'darkslategray4', '7' = 'deeppink4',
                   '9' = 'firebrick3')
# legend_colors = c('rho_2 = 20' = 'springgreen4', 'rho_2 = 25' = 'violetred4', 'rho_2 = 30' = 'slateblue4')
figure2a <- ggplot()
for (ri in 1:length(rho_2_check)) {

  # signal_name = paste0('rho2 = ', rho_2_check[ri])

  df_subset <- plot_df[plot_df$rho_2 == rho_2_check[ri], ]
  figure2a <- figure2a +
    geom_line(aes(x = eps, y = ci_width, color = rho_1_minus_rho_2), size = 1.0, alpha = 0.7, data = df_subset)
}
figure2a <- figure2a +
  xlab(TeX('$\\epsilon$')) + ylab('Average 90% CI width') +
  scale_color_manual(breaks = c('1', '3', '5', '7', '9'), values = legend_colors) +
  # ggtitle(TeX('Average confidence interval width as a function of $\\epsilon$')) +
  labs(color = TeX('$\\rho_1 - \\rho_2$')) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
figure2a
ggsave('figures/conf_width_epsilon_poisson.pdf',
       plot = figure2a, device = 'pdf', width = 4.5, height = 3.5)

# Figure 2b - Tradeoff of signal regime, RAND index for clustering, separated by epsilon
# --------------------------------------------------------------------------------------
#ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_true_index, K_index, , , ])
# ci_width <- apply(ci_width, c(1, 2), mean)
rand_slice <- apply(rand_results[n_index, K_true_index, K_index, , , ], c(1, 2), mean)
eps_index_plotting <- matrix(rep(eps_check, length(signal_regimes)), ncol = length(eps_check), byrow = TRUE)
rho_2_plotting <- matrix(rep(rho_2_check, length(eps_check)), ncol = length(eps_check))
rho_1_minus_rho_2_plotting <- matrix(rep(rho_1_minus_rho_2, length(eps_check)), ncol = length(eps_check))
signal_index_plotting <- matrix(rep(1:length(signal_regimes), length(eps_check)), ncol = length(eps_check))

legend_colors <- c('1' = 'gold3', '3' = 'darkseagreen4',
                   '5' = 'darkslategray4', '7' = 'deeppink4',
                   '9' = 'firebrick3')

plot_df <- data.frame(signal_index = as.vector(signal_index_plotting),
                      rho_2 = as.character(as.vector(rho_2_plotting)),
                      rho_1_minus_rho_2 = as.character(as.vector(rho_1_minus_rho_2_plotting)),
                      rand_index = as.vector(rand_slice),
                      eps = as.vector(eps_index_plotting))
figure2b <- ggplot(plot_df) +
  geom_line(aes(x = eps, y = rand_index, color = rho_1_minus_rho_2), size = 0.9) +
  scale_color_manual(breaks = c('100', '101', '102', '103', '104'), values = legend_colors) +
  xlab(TeX('$\\epsilon')) + ylab('Adjusted RAND Index') +
  coord_fixed() +
  #ggtitle(TeX('Classification accuracy (Adjusted RAND index) as a function of $\\epsilon$')) +
  labs(color = TeX('$\\rho_1 - \\rho_2$')) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
figure2b
ggsave('figures/randindex_epsilon_poisson.pdf',
       plot = figure2b, device = 'pdf', width = 4.5, height = 3.5)

