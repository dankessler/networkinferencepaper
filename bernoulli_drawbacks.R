library('networkinference')
library('latex2exp')
library('tidyverse')

ggplot2::theme_set(theme_minimal())

# ===========================
# == Simulation parameters ==
# ===========================

n <- 150                                         # Network nodes
K <- 2                                          # Number of communities
num_sim_per_gamma <- 200                        # Repetitions of the simulation
gamma_check <- seq(0.00, 0.50, length.out = 25) # What gamma to check
alpha <- 0.05                                   # Confidence level

# Marginal mean matrix
M <- matrix(rep(0.5, n^2), nrow = n)

# The linear combination vector we test
u <- c(1, -1, -1, 1)
u <- u / sqrt(sum(u^2))

# ============================================
# == Check through multiple values of gamma ==
# ============================================

# Where results are stored
targets <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
targets_tilde <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_midpoint <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_lower <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_upper <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]

  for (rep in 1:num_sim_per_gamma) {

    # Draw a matrix A from E[A] = M
    A <- matrix(rbinom(n^2, size = 1, prob = as.vector(M)), nrow = n)

    # Fission A
    W <- matrix(rbinom(n^2, size = 1, prob = gamma), nrow = n)
    A_tr <- A * (1 - W) + (1 - A) * W

    # ----------------
    # -- Clustering --
    # ----------------

    # Clustering on all of A
    z_hat_full_initial <- nett::spec_clust(A, K = K)
    z_hat_full <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat_full[, i] <- 1 * (z_hat_full_initial == i)
    }
    n_hat_full <- apply(z_hat_full, 2, sum)

    # Clustering on A_tr
    z_hat_initial <- nett::spec_clust(A_tr, K = K)
    z_hat <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat[, i] <- 1 * (z_hat_initial == i)
    }
    n_hat <- apply(z_hat, 2, sum)

    # ------------------
    # -- Do inference --
    # ------------------

    # Calculate our estimator and its target
    Cmask <- (gamma / (1-gamma))^(2*A_tr - 1)
    Tmat <- M / (M + (1-M) * Cmask)
    NN_inv <- diag(1 / diag(t(z_hat) %*% z_hat))

    Bhat <- NN_inv %*% t(z_hat) %*% A %*% z_hat %*% NN_inv
    B <- NN_inv %*% t(z_hat) %*% Tmat %*% z_hat %*% NN_inv
    Btilde <- NN_inv %*% t(z_hat) %*% M %*% z_hat %*% NN_inv

    # B_transform <- B / (B + (1 - B) * Cmask^(-1))
    # Bhat_transform <- Bhat / (Bhat + (1 - Bhat) * Cmask^(-1))

    Delta_hat <- matrix(0, nrow = K, ncol = K)

    for (k in 1:K) {
      for (l in 1:K) {
        k_nodes <- which(z_hat_initial == k)
        l_nodes <- which(z_hat_initial == l)

        # Variance
        Delta_hat[k, l] <- (Bhat[k, l] * (1 - Bhat[k, l])) / (length(k_nodes) * length(l_nodes))
      }
    }

    # Build a confidence interval
    alpha <- 0.10
    target <- t(u) %*% as.vector(B)
    target_marginal <- t(u) %*% as.vector(Btilde)
    estimate <- t(u) %*% as.vector(Bhat)

    estimate_var <- t(u) %*% diag(as.vector(Delta_hat)) %*% u
    margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)

    # Record the confidence interval
    targets[gamma_index, rep] <- target
    targets_tilde[gamma_index, rep] <- target_marginal
    ci_midpoint[gamma_index, rep] <- estimate
    ci_lower[gamma_index, rep] <- estimate - margin_of_error
    ci_upper[gamma_index, rep] <- estimate + margin_of_error
  }
}

# ====================
# == Create figures ==
# ====================

# -------------------------------------------------
# -- Confidence intervals as a function of gamma --
# -------------------------------------------------

avg_ci_midpoint <- apply(ci_midpoint, 1, mean)
avg_ci_lower <- apply(ci_lower, 1, mean)
avg_ci_upper <- apply(ci_upper, 1, mean)
avg_targets <- apply(targets, 1, mean)
avg_targets_tilde <- apply(targets_tilde, 1, mean)

plot_df <- data.frame(gamma = gamma_check,
                      avg_ci_midpoint = avg_ci_midpoint,
                      avg_ci_lower = avg_ci_lower,
                      avg_ci_upper = avg_ci_upper,
                      avg_targets = avg_targets,
                      avg_targets_tilde = avg_targets_tilde)

A_CI_gamma <- ggplot(plot_df) +
  geom_line(aes(x = gamma, y = avg_ci_lower, color = 'bounds'), linewidth = 1.1) +
  geom_line(aes(x = gamma, y = avg_ci_midpoint, color = 'estimate'), linewidth = 1.1) +
  geom_line(aes(x = gamma, y = avg_ci_upper, color = 'bounds'), linewidth = 1.1) +
  geom_line(aes(x = gamma, y = avg_targets, color = 'targets'), linewidth = 1.1) +
  geom_line(aes(x = gamma, y = avg_targets_tilde, color = 'targets_tilde'), linewidth = 1.1) +
  # geom_line(aes(x = gamma, y = 0), linetype = 'dashed', alpha = 0.6, linewidth = 0.8) +
  xlab(TeX('$\\gamma$')) + ylab('') + # ylab(TeX('$\\theta(A^{(tr)}_\\gamma)$')) +
  labs(color = 'Legend') +
  scale_color_manual(
    values = c('bounds' = 'darkslategray4', 'estimate' = 'firebrick3', 'targets' = 'darkorange', 'targets_tilde' = 'purple'),
    labels = c(TeX('CI lower/upper bounds'),
               TeX('$\\hat{\\theta}(A, A^{(tr)})$'),
               TeX('$\\theta(A^{(tr)})$'),
               TeX('$\\varphi(A^{(tr)})$'))
  ) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
A_CI_gamma
ggsave('figures/bernoulli_drawbacks/fig1.pdf', plot = A_CI_gamma,
       device = 'pdf', width = 6, height = 4.5)
