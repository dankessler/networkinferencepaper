# Libraries
library('ggplot2')
library('tidyverse')
library('nett')
library('mclust')
library('latex2exp')
library('patchwork')

ggplot2::theme_set(theme_minimal())

set.seed(1)

# ================
# == Simulation ==
# ================

num_sim <- 2000 # TODO: Raise this number
n <- 100
K <- 2
# gamma_check <- c(0.001, seq(0.05, 0.5, length.out = 10))
gamma_check <- c(0.001, 0.10, 0.20, 0.30, 0.40, 0.50)

# What contrast to test
u <- c(1, 0, 0, 0)
u <- u / sqrt(sum(u^2))

# Results of simulation
Bkl_results <- array(numeric(), dim = c(3, length(gamma_check), num_sim))
Vkl_results <- array(numeric(), dim = c(3, length(gamma_check), num_sim))
Phikl_results <- array(numeric(), dim = c(3, length(gamma_check), num_sim))

# Simulate
for (gamma_index in 1:length(gamma_check)) {
  cat(paste0('Gamma ', gamma_index, '/', length(gamma_check), '\n'))
  gamma <- gamma_check[gamma_index]

  # Settings 1, 2, and 3 for M
  M1 <- matrix(0.5, nrow = n, ncol = n)
  M2 <- rbind(cbind(matrix(0.6, nrow = n / 2, ncol = n / 2),
                    matrix(0.4, nrow = n / 2, ncol = n / 2)),
              cbind(matrix(0.4, nrow = n / 2, ncol = n / 2),
                    matrix(0.6, nrow = n / 2, ncol = n / 2)))
  # M3 <- matrix(sample(c(0.2, 0.8), n^2, rep = TRUE), nrow = n, ncol = n)
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

    # Repeatedly draw under this setting
    for (rep in 1:num_sim) {
      A <- matrix(rbinom(n^2, 1, as.vector(M)), nrow = n, ncol = n)
      W <- matrix(rbinom(n^2, 1, as.vector(rep(gamma, n^2))), nrow = n, ncol = n)
      A_tr <- A * (1-W) + (1-A) * W

      # Cluster based upon thinned matrix
      z_hat_initial <- nett::spec_clust(A_tr, K = K)
      z_hat <- matrix(rep(NA, n*K), nrow = n)
      for (i in 1:K) {
        z_hat[, i] <- 1 * (z_hat_initial == i)
      }
      n_hat <- apply(z_hat, 2, sum)

      # ---------------------
      # -- Useful matrices --
      # ---------------------

      NN_inv <- diag(1 / diag(t(z_hat) %*% z_hat))
      comm_pair_sample_size <- t(z_hat) %*% matrix(1, nrow = n, ncol = n) %*% z_hat
      comm_pair_sample_size_ones <- t(z_hat) %*% A_tr %*% z_hat
      comm_pair_sample_size_zeros <- t(z_hat) %*% (1 - A_tr) %*% z_hat

      # Weighting of ones and zeros
      ones_weighting <- comm_pair_sample_size_ones / comm_pair_sample_size
      zeros_weighting <- comm_pair_sample_size_zeros / comm_pair_sample_size

      # ---------------------------
      # -- Calculate all targets --
      # ---------------------------

      # B_{kl}
      Bkl <- t(u) %*% as.vector(NN_inv %*% t(z_hat) %*% M %*% z_hat %*% NN_inv)

      # V_{kl}
      Cmask <- (gamma / (1-gamma))^(2*A_tr - 1)
      Tmat <- M / (M + (1-M)*Cmask)
      c0 <- log(gamma / (1-gamma))
      c1 <- log((1-gamma) / gamma)

      Vkl <- t(u) %*% as.vector(NN_inv %*% t(z_hat) %*% Tmat %*% z_hat %*% NN_inv)

      # Phi_{kl}
      Lambda0_star <- (t(z_hat) %*% (Tmat * (1 - A_tr)) %*% z_hat) /
        comm_pair_sample_size_zeros
      Lambda1_star <- (t(z_hat) %*% (Tmat * A_tr) %*% z_hat) /
        comm_pair_sample_size_ones

      Phi_matrix <- zeros_weighting * (Lambda0_star / (Lambda0_star + (1 - Lambda0_star) * exp(c0))) +
        ones_weighting * (Lambda1_star / (Lambda1_star + (1 - Lambda1_star) * exp(c1)))

      Phikl <- t(u) %*% as.vector(Phi_matrix)

      # Store results
      Bkl_results[setting, gamma_index, rep] <- Bkl
      Vkl_results[setting, gamma_index, rep] <- Vkl
      Phikl_results[setting, gamma_index, rep] <- Phikl
    }
  }
}

# Average across the trials
avg_Bkl <- apply(Bkl_results, c(1, 2), mean)
avg_Vkl <- apply(Vkl_results, c(1, 2), mean)
avg_Phikl <- apply(Phikl_results, c(1, 2), mean)

# Legend labels
legend_values <- c('vkl_bkl' = 'deepskyblue3',
                   'phikl_bkl' = 'firebrick2')
legend_labels <- c('vkl_bkl' = TeX('$|V_{1,1} - B_{1,1}|$'),
                   'phikl_bkl' = TeX('$|\\Phi_{1,1} - B_{1,1}|$'))
shared_color_scale <- scale_color_manual(values = legend_values,
                                         labels = legend_labels,
                                         breaks = c('vkl_bkl', 'phikl_bkl'))

y_limits <- c(0, 0.015)

# Create a plot for each of the settings
settings_plots <- list()
for (setting in 1:3) {

  # Stack into appropriate data frame
  plot_df <- data.frame(gamma = gamma_check,
                        vkl_bkl = abs(avg_Bkl[setting, ] - avg_Vkl[setting, ]),
                        phikl_bkl = abs(avg_Bkl[setting, ] - avg_Phikl[setting, ]))
  plot_df <- plot_df %>%
    pivot_longer(c('vkl_bkl', 'phikl_bkl'))

  # Create the plot
  settings_plots[[setting]] <- ggplot(plot_df) +
    geom_line(aes(x = gamma, y = value, color = name), size = 0.8) +
    xlab(TeX('$\\gamma$')) + ylab('') +
    shared_color_scale +
    labs(color = 'Legend') +
    ylim(y_limits) +
    ggtitle(paste0('Setting ', setting), ) +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 10, angle = 0),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 11),
          # legend.position = 'none',
          legend.text = element_text(size = 10))
}

# With legend
row_plots <- (settings_plots[[1]] + settings_plots[[2]] + settings_plots[[3]]) +
  plot_layout(ncol = 3, guides = "collect")
row_plots


# Save plot
width <- 20
height <- width * (3.40 / 10)

ggsave('figures/bernoulli_fission_logic.pdf', plot = row_plots,
       device = 'pdf', width = width, height = height,
       units = 'cm')









# # Just plot each directly (no differences)

# Uncomment and this will run fine

# # Create a plot for each of the settings
# settings_plots <- list()
# for (setting in 1:3) {
#
#   # Stack into appropriate data frame
#   plot_df <- data.frame(gamma = gamma_check,
#                         vkl = avg_Vkl[setting, ],
#                         bkl = avg_Bkl[setting, ],
#                         phikl = avg_Phikl[setting, ])
#   plot_df <- plot_df %>%
#     pivot_longer(c('vkl', 'bkl', 'phikl'))
#
#   # Create the plot
#   settings_plots[[setting]] <- ggplot(plot_df) +
#     geom_line(aes(x = gamma, y = value, color = name)) +
#     xlab(TeX('$\\gamma$')) + ylab('') +
#     theme(aspect.ratio = 1,
#           axis.text.x = element_text(size = 16),
#           axis.text.y = element_text(size = 16),
#           axis.title.x = element_text(size = 25),
#           axis.title.y = element_text(size = 16),
#           legend.title = element_text(size = 18),
#           legend.text = element_text(size = 15))
# }
#
# settings_plots[[1]] + settings_plots[[2]] + settings_plots[[3]]
#
