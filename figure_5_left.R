# ==================================================== #
# == R code to replicate the left panel of Figure 5 == #
# == Author: Ethan Ancell                           == #
# ==================================================== #

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

n_check <- c(200) # Total network size
K_check <- c(2, 5, 10) # Guessed number of communities
K_true_check <- c(5)
eps_check <- c(0.5)
signal_regimes <- list(c(30, 27))

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

# ==============
# == Simulate ==
# ==============

cat("Simulating for Figure 5 - left panel\n")

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
            A <- matrix(rpois(n = n^2, lambda = as.vector(M)), nrow = n)

            # Split A
            A_split <- networkinference::split_network(A, distribution = "poisson",
                                                       epsilon = eps)
            Atr <- A_split$Atr
            Ate <- A_split$Ate

            # Clustering for naive and thinning
            z_hat <- nett::spec_clust(Atr, K = K)
            z_hat_naive <- nett::spec_clust(A, K = K)

            # Estimator, estimator variance, and target for the proposed approach
            thinning_infer <-
              networkinference::infer_network(Ate = Ate, u = u,
                                              communities = z_hat,
                                              distribution = "poisson",
                                              K = K, epsilon = eps)
            estimate <- thinning_infer$estimate
            estimate_var <- thinning_infer$estimate_variance
            target <- networkinference::check_target_of_inference(M = M, u = u,
                                                                  communities = z_hat,
                                                                  K = K)

            # Estimator, estimator variance, and target for the naive approach
            naive_infer <-
              networkinference::infer_network(Ate = A, u = u,
                                              communities = z_hat_naive,
                                              distribution = "poisson",
                                              K = K, epsilon = 0)
            estimate_naive <- naive_infer$estimate
            estimate_var_naive <- naive_infer$estimate_variance
            target_naive <- networkinference::check_target_of_inference(M = M, u = u,
                                                                        communities = z_hat_naive,
                                                                        K = K)

            # Save results
            # ------------
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

# ================================================ #
# == Save/load results (saves computation time) == #
# ================================================ #

# saveRDS(estimates_thinning, file = "saved_simulation_data/figure_5_left_estimates_thinning.RDS")
# saveRDS(targets_thinning, file = "saved_simulation_data/figure_5_left_targets_thinning.RDS")
# saveRDS(variances_thinning, file = "saved_simulation_data/figure_5_left_variances_thinning.RDS")
# saveRDS(estimates_naive, file = "saved_simulation_data/figure_5_left_estimates_naive.RDS")
# saveRDS(targets_naive, file = "saved_simulation_data/figure_5_left_targets_naive.RDS")
# saveRDS(variances_naive, file = "saved_simulation_data/figure_5_left_variances_naive.RDS")

estimates_thinning <- readRDS("saved_simulation_data/figure_5_left_estimates_thinning.RDS")
targets_thinning <- readRDS("saved_simulation_data/figure_5_left_targets_thinning.RDS")
variances_thinning <- readRDS("saved_simulation_data/figure_5_left_variances_thinning.RDS")
estimates_naive <- readRDS("saved_simulation_data/figure_5_left_estimates_naive.RDS")
targets_naive <- readRDS("saved_simulation_data/figure_5_left_targets_naive.RDS")
variances_naive <- readRDS("saved_simulation_data/figure_5_left_variances_naive.RDS")

# ========== #
# == Plot == #
# ========== #

checked_coverages <- seq(0, 1.0, length.out = 100)
n_coverages <- length(checked_coverages)

n_index <- 1
K_true_index <- 1
K_true <- K_true_check[K_true_index]
signal_index <- 1

# Results
figure1_coverage_thinning <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))
figure1_coverage_naive <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))

# Calculate coverages
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

# Stack all this information in a plotting data frame
plot_df <- data.frame()
for (K_index in 1:length(K_check)) {
  for (eps_index in 1:length(eps_check)) {
    # Proposed approach
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_thinning[K_index, eps_index, ],
                                method = rep('Proposed', n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages)))

    # Naive approach
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_naive[K_index, eps_index, ],
                                method = rep('Naive', n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages)))
  }
}

# Combine to create plot
legend_colors <- c('2' = 'cadetblue3', '5' = 'royalblue3', '10' = 'blue4')
figure <- ggplot()
for (K_index in 1:length(K_check)) {
  K <- K_check[K_index]

  # Subset data
  df_subset_thin <- plot_df %>% filter(method == 'Proposed', K == K)
  df_subset_naive <- plot_df %>% filter(method == 'Naive', K == K)

  # Add it onto the plot
  figure <- figure +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Proposed'), linewidth = 0.9, data = df_subset_thin)
  figure <- figure +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Naive'), linewidth = 0.7, data = df_subset_naive)
}
figure <- figure +
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

figure
ggsave('figures/nominal_empirical_coverage_poisson.pdf',
       plot = figure, device = 'pdf', width = 4.5, height = 3.5)
