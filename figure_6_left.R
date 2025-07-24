# ==================================================== #
# == R code to replicate the left panel of Figure 6 == #
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
gamma_check <- c(0.25)
signal_regimes <- list(c(0.75, 0.5))

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

# ==============
# == Simulate ==
# ==============

cat("Simulating for Figure 6 - left panel\n")

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
                                              K = K, gamma = gamma,
                                              Atr = Atr)
            estimate <- fission_infer$estimate
            estimate_var <- fission_infer$estimate_variance
            target_theta <- networkinference::check_target_of_inference(M = M, u = u,
                                                                  communities = z_hat,
                                                                  K = K)
            target_xi <- networkinference::check_target_of_inference(M = M, u = u,
                                                                     communities = z_hat,
                                                                     K = K,
                                                                     bernoulli_target = TRUE,
                                                                     gamma = gamma,
                                                                     Atr = Atr)

            # Estimator, estimator variance, and target for the naive approach
            # (Note that this is more complicated of code because my R package
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

            # Save results
            # ------------
            xi_est_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              estimate
            xi_target_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              target_xi
            theta_target_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              target_theta
            xi_variances_thinning[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              estimate_var

            theta_est_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              theta_est_naive_
            theta_target_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              theta_naive_
            theta_variances_naive[n_index, K_true_index, K_index, signal_regime_index, gamma_index, rep] <-
              theta_var_naive_
          }
        }
      }
    }
  }
}

# ================================================ #
# == Save/load results (saves computation time) == #
# ================================================ #

# saveRDS(xi_est_thinning, file = "saved_simulation_data/figure_6_left_xi_est_thinning.RDS")
# saveRDS(xi_target_thinning, file = "saved_simulation_data/figure_6_left_xi_target_thinning.RDS")
# saveRDS(theta_target_thinning, file = "saved_simulation_data/figure_6_left_theta_target_thinning.RDS")
# saveRDS(xi_variances_thinning, file = "saved_simulation_data/figure_6_left_xi_variances_thinning.RDS")
# saveRDS(theta_est_naive, file = "saved_simulation_data/figure_6_left_theta_est_naive.RDS")
# saveRDS(theta_target_naive, file = "saved_simulation_data/figure_6_left_theta_target_naive.RDS")
# saveRDS(theta_variances_naive, file = "saved_simulation_data/figure_6_left_theta_variances_naive.RDS")

xi_est_thinning <- readRDS("saved_simulation_data/figure_6_left_xi_est_thinning.RDS")
xi_target_thinning <- readRDS("saved_simulation_data/figure_6_left_xi_target_thinning.RDS")
theta_target_thinning <- readRDS("saved_simulation_data/figure_6_left_theta_target_thinning.RDS")
xi_variances_thinning <- readRDS("saved_simulation_data/figure_6_left_xi_variances_thinning.RDS")
theta_est_naive <- readRDS("saved_simulation_data/figure_6_left_theta_est_naive.RDS")
theta_target_naive <- readRDS("saved_simulation_data/figure_6_left_theta_target_naive.RDS")
theta_variances_naive <- readRDS("saved_simulation_data/figure_6_left_theta_variances_naive.RDS")

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
figure_coverage_thinning <- array(numeric(), c(length(K_check), length(gamma_check), length(checked_coverages)))
figure_coverage_naive <- array(numeric(), c(length(K_check), length(gamma_check), length(checked_coverages)))

for (K_index in 1:length(K_check)) {
  for (gamma_index in 1:length(gamma_check)) {
    for (coverage_index in 1:length(checked_coverages)) {

      coverage <- checked_coverages[coverage_index]
      z_quantile <- qnorm(1 - (1 - coverage) / 2)

      # Extract relevant slices from the results array
      targets_thinning_slice <- theta_target_thinning[n_index, K_true_index, K_index, signal_index, gamma_index, ]
      theta_target_naive_slice <- theta_target_naive[n_index, K_true_index, K_index, signal_index, gamma_index, ]

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
      figure_coverage_thinning[K_index, gamma_index, coverage_index] <-
        mean((lower_bounds_thinning <= targets_thinning_slice) & (targets_thinning_slice <= upper_bounds_thinning))
      figure_coverage_naive[K_index, gamma_index, coverage_index] <-
        mean((lower_bounds_naive <= theta_target_naive_slice) & (theta_target_naive_slice <= upper_bounds_naive))
    }
  }
}

# Stack all this information in a plotting data frame
plot_df <- data.frame()
for (K_index in 1:length(K_check)) {
  for (gamma_index in 1:length(gamma_check)) {
    # Thinning
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure_coverage_thinning[K_index, gamma_index, ],
                                method = rep('Proposed', n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                gamma = rep(gamma_check[gamma_index], n_coverages)))

    # Naive
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure_coverage_naive[K_index, gamma_index, ],
                                method = rep('Naive', n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                gamma = rep(gamma_check[gamma_index], n_coverages)))
  }
}

# Create plot
legend_colors <- c('2' = 'cadetblue3', '5' = 'royalblue3', '10' = 'blue4')
figure <- ggplot()
for (K_index in 1:length(K_check)) {
  K <- K_check[K_index]

  # Subset data
  df_subset_thin <- plot_df %>% filter(method == 'Proposed', K == K)
  df_subset_naive <- plot_df %>% filter(method == 'Naive', K == K)

  # Add it onto the plot
  figure <- figure +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Proposed'), size = 0.9, data = df_subset_thin)
  figure <- figure +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Naive'), size = 0.7, data = df_subset_naive)
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
ggsave('figures/nominal_empirical_coverage_bernoulli.pdf',
       plot = figure, device = 'pdf', width = 4.5, height = 3.5)
