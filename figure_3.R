# ================================== #
# == R code to replicate Figure 3 == #
# == Author: Ethan Ancell         == #
# ================================== #

library("ggplot2")
library("tidyverse")
library("nett")
library("mclust")
library("latex2exp")
library("patchwork")
library("networkinference")

ggplot2::theme_set(theme_minimal())
set.seed(1)

# ================
# == Simulation ==
# ================

num_sim <- 5000
n <- 100
K <- 2
gamma_check <- c(0.001, 0.10, 0.20, 0.30, 0.40, 0.50)

# Linear combination to check
u <- c(1, 0, 0, 0)

# Results of simulation
B_results <- array(numeric(), dim = c(3, length(gamma_check), num_sim))
V_results <- array(numeric(), dim = c(3, length(gamma_check), num_sim))
Phi_results <- array(numeric(), dim = c(3, length(gamma_check), num_sim))

# Simulate
for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]
  cat(paste0('Gamma = ', gamma, " (", gamma_index, ' out of ', length(gamma_check), ')\n'))

  # Settings 1, 2, and 3 for M
  M1 <- matrix(0.5, nrow = n, ncol = n)
  M2 <- rbind(cbind(matrix(0.6, nrow = n / 2, ncol = n / 2),
                    matrix(0.4, nrow = n / 2, ncol = n / 2)),
              cbind(matrix(0.4, nrow = n / 2, ncol = n / 2),
                    matrix(0.6, nrow = n / 2, ncol = n / 2)))
  M3 <- matrix(runif(n^2), nrow = n, ncol = n)

  settings <- c(1, 2, 3)
  for (setting in settings) {
    if (setting == 1) {
      M <- M1
    } else if (setting == 2) {
      M <- M2
    } else {
      M <- M3
    }

    # Repeatedly draw under this setting a total of num_sim times
    for (rep in 1:num_sim) {
      # Draw and fission matrix
      A <- matrix(rbinom(n^2, 1, as.vector(M)), nrow = n, ncol = n)

      split_A <- networkinference::split_matrix(A, distribution = "bernoulli",
                                                gamma = gamma)
      Atr <- split_A$Atr

      # Cluster with Atr
      z_hat <- nett::spec_clust(Atr, K = K)

      # Calculate targets
      B <- networkinference::check_target_of_inference(M = M, u = u,
                                                       communities = z_hat,
                                                       K = K)
      Phi <- networkinference::check_target_of_inference(M = M, u = u,
                                                         communities = z_hat,
                                                         K = K,
                                                         bernoulli_target = TRUE,
                                                         gamma = gamma,
                                                         Atr = Atr)

      # Manually calculate V since my R package does not provide functionality for this
      z_hat_matrix <- matrix(rep(NA, n*K), nrow = n)
      for (i in 1:K) {
        z_hat_matrix[, i] <- 1 * (z_hat == i)
      }
      n_hat <- apply(z_hat_matrix, 2, sum)
      NN_inv <- diag(1 / diag(t(z_hat_matrix) %*% z_hat_matrix))
      Cmask <- (gamma / (1-gamma))^(2*Atr - 1)
      Tmat <- M / (M + (1-M)*Cmask)

      V <- t(u) %*% as.vector(NN_inv %*% t(z_hat_matrix) %*% Tmat %*% z_hat_matrix %*% NN_inv)

      # Store results
      B_results[setting, gamma_index, rep] <- B
      V_results[setting, gamma_index, rep] <- V
      Phi_results[setting, gamma_index, rep] <- Phi
    }
  }
}

# ================================================ #
# == Save/load results (saves computation time) == #
# ================================================ #

# saveRDS(B_results, file = "saved_simulation_data/figure_3_B_results.RDS")
# saveRDS(V_results, file = "saved_simulation_data/figure_3_V_results.RDS")
# saveRDS(Phi_results, file = "saved_simulation_data/figure_3_Phi_results.RDS")

B_results <- readRDS("saved_simulation_data/figure_3_B_results.RDS")
V_results <- readRDS("saved_simulation_data/figure_3_V_results.RDS")
Phi_results <- readRDS("saved_simulation_data/figure_3_Phi_results.RDS")

# ========== #
# == Plot == #
# ========== #

# Take average of the differences between simulation settings
avg_abs_V_minus_B <- apply(abs(V_results - B_results), c(1, 2), mean)
avg_abs_Phi_minus_B <- apply(abs(Phi_results - B_results), c(1, 2), mean)

# Legend labels
legend_values <- c('V_minus_B' = 'deepskyblue3',
                   'Phi_minus_B' = 'firebrick2')
legend_labels <- c('V_minus_B' = TeX('$|V_{11}(A^{(tr)}) - B_{11}(A^{(tr)})|$'),
                   'Phi_minus_B' = TeX('$|\\Phi_{11}(A^{(tr)}) - B_{11}(A^{(tr)})|$'))
shared_color_scale <- scale_color_manual(values = legend_values,
                                         labels = legend_labels,
                                         breaks = c('V_minus_B', 'Phi_minus_B'))

y_limits <- c(0, 0.025)

# Create a plot for each of the settings
settings_plots <- list()
for (setting in 1:3) {

  # Stack into appropriate data frame
  plot_df <- data.frame(gamma = gamma_check,
                        V_minus_B = avg_abs_V_minus_B[setting, ],
                        Phi_minus_B = avg_abs_Phi_minus_B[setting, ])
  plot_df <- plot_df %>%
    pivot_longer(c('V_minus_B', 'Phi_minus_B'))

  # Create the plot
  settings_plots[[setting]] <- ggplot(plot_df) +
    geom_line(aes(x = gamma, y = value, color = name), linewidth = 0.8) +
    xlab(TeX('$\\gamma$')) + ylab('') +
    shared_color_scale +
    labs(color = 'Legend') +
    ylim(y_limits) +
    ggtitle(paste0('Setting ', setting), ) +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 7.8, angle = 0),
          axis.text.y = element_text(size = 7.8),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10))
}

# With legend
row_plots <- (settings_plots[[1]] + settings_plots[[2]] + settings_plots[[3]]) +
  plot_layout(ncol = 3, guides = "collect")
row_plots

# Save plot
width <- 20
height <- width * (2.9 / 10)

ggsave('figures/bernoulli_fission_logic.pdf', plot = row_plots,
       device = 'pdf', width = width, height = height,
       units = 'cm')
