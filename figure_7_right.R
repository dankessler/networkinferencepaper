# ===================================================== #
# == R code to replicate the right panel of Figure 7 == #
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
K_check <- c(5) # Guessed number of communities
gamma_check <- c(0.001, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20,
                 0.30, 0.40, 0.499)

# ====================== #
# == Simulation setup == #
# ====================== #

xi_est_thinning <- array(0, dim = c(length(n_check),
                                    length(K_check),
                                    length(gamma_check),
                                    num_sim))
xi_target_thinning <- array(0, dim = c(length(n_check),
                                       length(K_check),
                                       length(gamma_check),
                                       num_sim))
theta_target_thinning <- array(0, dim = c(length(n_check),
                                          length(K_check),
                                          length(gamma_check),
                                          num_sim))
xi_variances_thinning <- array(0, dim = c(length(n_check),
                                          length(K_check),
                                          length(gamma_check),
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

    # Different values of gamma in fission
    for (gamma_index in 1:length(gamma_check)) {
      gamma <- gamma_check[gamma_index]
      cat(paste0("\t\t(gamma = ", gamma, ")\n"))

      # -------------------------------------------------
      # -- Actual simulation part (non just for-loops) --
      # -------------------------------------------------

      # No such thing as real communities
      M <- matrix(runif(n^2), nrow = n, ncol = n)

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

        # Save results
        # ------------
        xi_est_thinning[n_index, K_index, gamma_index, rep] <-
          estimate
        xi_target_thinning[n_index, K_index, gamma_index, rep] <-
          target_xi
        theta_target_thinning[n_index, K_index, gamma_index, rep] <-
          target_theta
        xi_variances_thinning[n_index, K_index, gamma_index, rep] <-
          estimate_var
      }
    }
  }
}

# ================================================ #
# == Save/load results (saves computation time) == #
# ================================================ #

# saveRDS(xi_est_thinning, file = "saved_simulation_data/figure_7_right_xi_est_thinning.RDS")
# saveRDS(xi_target_thinning, file = "saved_simulation_data/figure_7_right_xi_target_thinning.RDS")
# saveRDS(theta_target_thinning, file = "saved_simulation_data/figure_7_right_theta_target_thinning.RDS")
# saveRDS(xi_variances_thinning, file = "saved_simulation_data/figure_7_right_xi_variances_thinning.RDS")

xi_est_thinning <- readRDS("saved_simulation_data/figure_7_right_xi_est_thinning.RDS")
xi_target_thinning <- readRDS("saved_simulation_data/figure_7_right_xi_target_thinning.RDS")
theta_target_thinning <- readRDS("saved_simulation_data/figure_7_right_theta_target_thinning.RDS")
xi_variances_thinning <- readRDS("saved_simulation_data/figure_7_right_xi_variances_thinning.RDS")

# ========== #
# == Plot == #
# ========== #

alpha <- 0.10
ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(xi_variances_thinning[n_index, K_index, , ])
ci_width <- apply(ci_width, 1, mean)

plot_df <- data.frame(ci_width = ci_width,
                      gamma = gamma_check,
                      one_minus_gamma = 1 - gamma_check)

figure <- ggplot(data = plot_df) +
  geom_line(aes(x = one_minus_gamma, y = ci_width), linewidth = 1.0, alpha = 0.7)

figure <- figure +
  xlab(TeX('$1-\\gamma$')) + ylab('Average 90% CI width') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16))
figure
ggsave('figures/conf_width_gamma_bernoulli_maximal_heterogeneity.pdf',
       plot = figure, device = 'pdf', width = 4.5, height = 3.5)
