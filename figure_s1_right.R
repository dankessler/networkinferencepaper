# ====================================================== #
# == R code to replicate the right panel of Figure S1 == #
# == Author: Ethan Ancell                             == #
# ====================================================== #

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
eps_check <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
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

        # Split A using the split_network() function
        A_split <- networkinference::split_network(A, distribution = 'gaussian',
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

# saveRDS(estimates_thinning, file = "saved_simulation_data/figure_s1_right_estimates_thinning.RDS")
# saveRDS(targets_thinning, file = "saved_simulation_data/figure_s1_right_targets_thinning.RDS")
# saveRDS(variances_thinning, file = "saved_simulation_data/figure_s1_right_variances_thinning.RDS")
# saveRDS(estimates_naive, file = "saved_simulation_data/figure_s1_right_estimates_naive.RDS")
# saveRDS(targets_naive, file = "saved_simulation_data/figure_s1_right_targets_naive.RDS")
# saveRDS(variances_naive, file = "saved_simulation_data/figure_s1_right_variances_naive.RDS")

estimates_thinning <- readRDS("saved_simulation_data/figure_s1_right_estimates_thinning.RDS")
targets_thinning <- readRDS("saved_simulation_data/figure_s1_right_targets_thinning.RDS")
variances_thinning <- readRDS("saved_simulation_data/figure_s1_right_variances_thinning.RDS")
estimates_naive <- readRDS("saved_simulation_data/figure_s1_right_estimates_naive.RDS")
targets_naive <- readRDS("saved_simulation_data/figure_s1_right_targets_naive.RDS")
variances_naive <- readRDS("saved_simulation_data/figure_s1_right_variances_naive.RDS")

# ========== #
# == Plot == #
# ========== #

n_index <- 1
K_index <- 1
alpha <- 0.10

ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_index, , ])
ci_width <- apply(ci_width, 1, mean)

plot_df <- data.frame(ci_width = ci_width,
                      eps = eps_check)

figure <- ggplot(data = plot_df) +
  geom_line(aes(x = eps, y = ci_width), linewidth = 1.0, alpha = 0.7) +
  xlab(TeX('$\\epsilon$')) + ylab('Average 90% CI width') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16))
figure
ggsave('figures/gaussian_heterogeneous_ciwidth.pdf',
       plot = figure, device = 'pdf', width = 4.5, height = 3.5)
