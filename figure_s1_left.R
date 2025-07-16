# ===================================================== #
# == R code to replicate the left panel of Figure S1 == #
# == Author: Ethan Ancell                            == #
# ===================================================== #

library("networkinference")
library("ggplot2")
library("tidyverse")
library("nett")
library("mclust")
library("latex2exp")

ggplot2::theme_set(theme_minimal())

# ========================= #
# == Simulation settings == #
# ========================= #

set.seed(1)
num_sim <- 5000 # Number of replications
n_check <- c(200) # Total network size
K_check <- c(2, 5, 10) # Guessed number of communities
eps_check <- c(0.50)
tau2 <- 25

# ====================== #
# == Simulation setup == #
# ====================== #

estimates_thinning <- array(0, dim = c(length(n_check),
                                       length(K_check),
                                       length(eps_check),
                                       num_sim))
targets_thinning <- array(0, dim = c(length(n_check),
                                     length(K_check),
                                     length(eps_check),
                                     num_sim))
variances_thinning <- array(0, dim = c(length(n_check),
                                       length(K_check),
                                       length(eps_check),
                                       num_sim))

estimates_naive <- array(0, dim = c(length(n_check),
                                    length(K_check),
                                    length(eps_check),
                                    num_sim))
targets_naive <- array(0, dim = c(length(n_check),
                                  length(K_check),
                                  length(eps_check),
                                  num_sim))
variances_naive <- array(0, dim = c(length(n_check),
                                    length(K_check),
                                    length(eps_check),
                                    num_sim))

# ==============
# == Simulate ==
# ==============

# Sample size
for (n_index in 1:length(n_check)) {
  n <- n_check[n_index]
  cat(paste0("(n = ", n, ")\n"))

  # Number of communities used in clustering
  for (K_index in 1:length(K_check)) {
    K <- K_check[K_index]
    cat(paste0("\t(K = ", K, ")\n"))

    # Different values of epsilon in thinning
    for (eps_index in 1:length(eps_check)) {
      eps <- eps_check[eps_index]
      cat(paste0('\t\t(epsilon = ', eps, ')\n'))

      # -------------------------------------------------
      # -- Actual simulation part (non just for-loops) --
      # -------------------------------------------------

      # No such thing as real communities
      M <- matrix(runif(n^2, min = 0, max = 20), nrow = n, ncol = n)

      for (rep in 1:num_sim) {

        # Linear combination vector
        u <- c(1, rep(0, K^2 - 1))

        # Draw A
        A <- matrix(rnorm(n = n^2, mean = as.vector(M), sd = sqrt(tau2)), nrow = n)

        # Split A using the split_matrix() function
        A_split <- networkinference::split_matrix(A, distribution = 'gaussian',
                                                  epsilon = eps, allow_self_loops = TRUE,
                                                  is_directed = TRUE, tau = sqrt(tau2))
        A_tr <- A_split$Atr
        A_te <- A_split$Ate

        # Cluster (thinning)
        z_hat_initial <- nett::spec_clust(A_tr, K = K)
        z_hat <- matrix(rep(NA, n*K), nrow = n)
        for (i in 1:K) {
          z_hat[, i] <- 1 * (z_hat_initial == i)
        }
        n_hat <- apply(z_hat, 2, sum)

        # Cluster (naive)
        z_hat_initial_naive <- nett::spec_clust(A, K = K)
        z_hat_naive <- matrix(rep(NA, n*K), nrow = n)
        for (i in 1:K) {
          z_hat_naive[, i] <- 1 * (z_hat_initial_naive == i)
        }
        n_hat_naive <- apply(z_hat_naive, 2, sum)

        # Proposed approach - Conduct inference
        inference_result <-
          networkinference::infer_network(Ate = A_te, u = u,
                                          communities = z_hat_initial,
                                          distribution = 'gaussian',
                                          epsilon = eps, K = K,
                                          tau = sqrt(tau2))
        estimate <- inference_result$estimate
        estimate_var <- inference_result$estimate_variance

        target <-
          networkinference::check_target_of_inference(M = M, u = u,
                                                      communities = z_hat_initial,
                                                      K = K)

        # Naive approach - Conduct inference
        naive_inference_result <-
          networkinference::infer_network(Ate = A, u = u,
                                          communities = z_hat_initial_naive,
                                          distribution = 'gaussian',
                                          epsilon = 0, K = K,
                                          tau = sqrt(tau2))
        estimate_naive <- naive_inference_result$estimate
        estimate_var_naive <- naive_inference_result$estimate_var

        target_naive <-
          networkinference::check_target_of_inference(M = M, u = u,
                                                      communities = z_hat_initial_naive,
                                                      K = K)

        # Save everything else
        # --------------------
        estimates_thinning[n_index, K_index, eps_index, rep] <-
          estimate
        targets_thinning[n_index, K_index, eps_index, rep] <-
          target
        variances_thinning[n_index, K_index, eps_index, rep] <-
          estimate_var
        estimates_naive[n_index, K_index, eps_index, rep] <-
          estimate_naive
        targets_naive[n_index, K_index, eps_index, rep] <-
          target_naive
        variances_naive[n_index, K_index, eps_index, rep] <-
          estimate_var_naive
      }
    }
  }
}

# ================================================ #
# == Save/load results (saves computation time) == #
# ================================================ #

saveRDS(estimates_thinning, file = "saved_simulation_data/figure_s1_left_estimates_thinning.RDS")
saveRDS(targets_thinning, file = "saved_simulation_data/figure_s1_left_targets_thinning.RDS")
saveRDS(variances_thinning, file = "saved_simulation_data/figure_s1_left_variances_thinning.RDS")
saveRDS(estimates_naive, file = "saved_simulation_data/figure_s1_left_estimates_naive.RDS")
saveRDS(targets_naive, file = "saved_simulation_data/figure_s1_left_targets_naive.RDS")
saveRDS(variances_naive, file = "saved_simulation_data/figure_s1_left_variances_naive.RDS")

estimates_thinning <- readRDS("saved_simulation_data/figure_s1_left_estimates_thinning.RDS")
targets_thinning <- readRDS("saved_simulation_data/figure_s1_left_targets_thinningRDS")
variances_thinning <- readRDS("saved_simulation_data/figure_s1_left_variances_thinning.RDS")
estimates_naive <- readRDS("saved_simulation_data/figure_s1_left_estimates_naive.RDS")
targets_naive <- readRDS("saved_simulation_data/figure_s1_left_targets_naive.RDS")
variances_naive <- readRDS("saved_simulation_data/figure_s1_left_variances_naive.RDS")

# ========== #
# == Plot == #
# ========== #

checked_coverages <- seq(0, 1.0, length.out = 100)
n_coverages <- length(checked_coverages)

n_index <- 1

# Results
figure_coverage_thinning <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))
figure_coverage_naive <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))

for (K_index in 1:length(K_check)) {

  for (eps_index in 1:length(eps_check)) {

    for (coverage_index in 1:length(checked_coverages)) {

      coverage <- checked_coverages[coverage_index]
      z_quantile <- qnorm(1 - (1 - coverage) / 2)

      # Extract relevant slices from the results array
      targets_thinning_slice <- targets_thinning[n_index, K_index, eps_index, ]
      theta_target_naive_slice <- targets_naive[n_index, K_index, eps_index, ]

      upper_bounds_thinning <-
        estimates_thinning[n_index, K_index, eps_index, ] +
        z_quantile * sqrt(variances_thinning[n_index, K_index, eps_index, ])
      lower_bounds_thinning <-
        estimates_thinning[n_index, K_index, eps_index, ] -
        z_quantile * sqrt(variances_thinning[n_index, K_index, eps_index, ])

      upper_bounds_naive <-
        estimates_naive[n_index, K_index, eps_index, ] +
        z_quantile * sqrt(variances_naive[n_index, K_index, eps_index, ])
      lower_bounds_naive <-
        estimates_naive[n_index, K_index, eps_index, ] -
        z_quantile * sqrt(variances_naive[n_index, K_index, eps_index, ])

      # Store coverage
      figure_coverage_thinning[K_index, eps_index, coverage_index] <-
        mean((lower_bounds_thinning <= targets_thinning_slice) & (targets_thinning_slice <= upper_bounds_thinning))
      figure_coverage_naive[K_index, eps_index, coverage_index] <-
        mean((lower_bounds_naive <= theta_target_naive_slice) & (theta_target_naive_slice <= upper_bounds_naive))

    }

  }

}

# Stack all this information in a data frame that can be used for plotting
plot_df <- data.frame()
for (K_index in 1:length(K_check)) {
  for (eps_index in 1:length(eps_check)) {

    # Thinning
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure_coverage_thinning[K_index, eps_index, ],
                                method = rep('Proposed', n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages),
                                eps_ch = as.character(rep(eps_check[eps_index], n_coverages))))

    # Naive
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure_coverage_naive[K_index, eps_index, ],
                                method = rep('Naive', n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages),
                                eps_ch = as.character(rep(eps_check[eps_index], n_coverages))))
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
ggsave('figures/gaussian_heterogeneous_coverage.pdf',
       plot = figure, device = 'pdf', width = 4.5, height = 3.5)
