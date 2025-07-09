# ============================================
# == Simulations for network thinning paper ==
# == (Binary edges)                         ==
# == Author: Ethan Ancell                   ==
# ============================================

# Libraries
library('ggplot2')
library('tidyverse')
library('nett')
library('mclust')
library('latex2exp')

ggplot2::theme_set(theme_minimal())

# ---------------------------------------------- #
# -- Average confidence interval width figure -- #
# ---------------------------------------------- #

set.seed(1)
num_sim <- 5000 # Number of replications

n_check <- c(200) # Total network size
K_check <- c(5) # Guessed number of communities
eps_check <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
tau2 <- 25

alpha <- 0.10
use_random_u <- FALSE

# Where results are stored
# ------------------------
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

# Coverage results of (1-alpha) confidence intervals
covered_thinning <- array(0, dim = c(length(n_check),
                                     length(K_check),
                                     length(eps_check)))
covered_naive <- array(0, dim = c(length(n_check),
                                  length(K_check),
                                  length(eps_check)))

# ==============
# == Simulate ==
# ==============

# Sample size
for (n_index in 1:length(n_check)) {
  n <- n_check[n_index]
  cat(paste0('(n = ', n, ')\n'))

  # Number of communities used in clustering
  for (K_index in 1:length(K_check)) {
    K <- K_check[K_index]
    cat(paste0('\t(K = ', K, ')\n'))

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

        # Sample u
        if (use_random_u) {
          # Sample the vector u from the unit K^2-sphere
          u <- rnorm(K^2)
          u <- u / sqrt(sum(u^2))
        } else {
          u <- c(1, rep(0, K^2 - 1))
        }

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

        # Construct a confidence interval for the target of inference
        margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)
        margin_of_error_naive <- qnorm(1 - alpha / 2) * sqrt(estimate_var_naive)

        # Check coverage (thinning and naive)
        covered_thinning[n_index, K_index, eps_index] <-
          covered_thinning[n_index, K_index, eps_index] +
          as.integer((estimate - margin_of_error <= target) & (target <= estimate + margin_of_error))

        covered_naive[n_index, K_index, eps_index] <-
          covered_naive[n_index, K_index, eps_index] +
          as.integer((estimate_naive - margin_of_error_naive <= target_naive) & (target_naive <= estimate_naive + margin_of_error_naive))

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

covered_thinning <- covered_thinning / num_sim
covered_naive <- covered_naive / num_sim

covered_thinning
covered_naive

# ----------------------------------------------
# - Save/load results (saves computation time) -
# ----------------------------------------------

# saveRDS(estimates_thinning, file = 'save_data/gaussian_heterogeneous_ciwidth_estimates_thinning.RDS')
# saveRDS(targets_thinning, file = 'save_data/gaussian_heterogeneous_ciwidth_targets_thinning.RDS')
# saveRDS(variances_thinning, file = 'save_data/gaussian_heterogeneous_ciwidth_variances_thinning.RDS')
# saveRDS(estimates_naive, file = 'save_data/gaussian_heterogeneous_ciwidth_estimates_naive.RDS')
# saveRDS(targets_naive, file = 'save_data/gaussian_heterogeneous_ciwidth_targets_naive.RDS')
# saveRDS(variances_naive, file = 'save_data/gaussian_heterogeneous_ciwidth_variances_naive.RDS')
# saveRDS(covered_thinning, file = 'save_data/gaussian_heterogeneous_ciwidth_covered_thinning.RDS')
# saveRDS(covered_naive, file = 'save_data/gaussian_heterogeneous_ciwidth_covered_naive.RDS')

estimates_thinning <- readRDS('save_data/gaussian_heterogeneous_ciwidth_estimates_thinning.RDS')
targets_thinning <- readRDS('save_data/gaussian_heterogeneous_ciwidth_targets_thinning.RDS')
variances_thinning <- readRDS('save_data/gaussian_heterogeneous_ciwidth_variances_thinning.RDS')
estimates_naive <- readRDS('save_data/gaussian_heterogeneous_ciwidth_estimates_naive.RDS')
targets_naive <- readRDS('save_data/gaussian_heterogeneous_ciwidth_targets_naive.RDS')
variances_naive <- readRDS('save_data/gaussian_heterogeneous_ciwidth_variances_naive.RDS')
covered_thinning <- readRDS('save_data/gaussian_heterogeneous_ciwidth_covered_thinning.RDS')
covered_naive <- readRDS('save_data/gaussian_heterogeneous_ciwidth_covered_naive.RDS')

# ============== #
# == Plotting == #
# ============== #

# ----------------------------- #
# - Confidence interval width - #
# ----------------------------- #

n_index <- 1
K_index <- 1
alpha <- 0.10

# Figure 2a - CI width as a function of signal, separate by eps
# -----------------------------------------------------------------
ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_index, , ])
ci_width <- apply(ci_width, 1, mean)

plot_df <- data.frame(ci_width = ci_width,
                      eps = eps_check)

figure2a <- ggplot(data = plot_df) +
  geom_line(aes(x = eps, y = ci_width), linewidth = 1.0, alpha = 0.7) +
  xlab(TeX('$\\epsilon$')) + ylab('Average 90% CI width') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16))
figure2a
ggsave('figures/gaussian_heterogeneous_ciwidth.pdf',
       plot = figure2a, device = 'pdf', width = 4.5, height = 3.5)
