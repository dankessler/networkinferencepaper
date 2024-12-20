library('networkinference')
library('latex2exp')
library('tidyverse')

ggplot2::theme_set(theme_minimal())

# =====================================================
# == Manually load the Zachary's karate club dataset ==
# =====================================================

n <- 34
K <- 2

zachary <- matrix(0, nrow = n, ncol = n)

# List of connections
connections <- list(
  c(2, 1), c(3, 1), c(3, 2), c(4, 1), c(4, 2), c(4, 3), c(5, 1), c(6, 1),
  c(7, 1), c(7, 5), c(7, 6), c(8, 1), c(8, 2), c(8, 3), c(8, 4), c(9, 1),
  c(9, 3), c(10, 3), c(11, 1), c(11, 5), c(11, 6), c(12, 1), c(13, 1),
  c(13, 4), c(14, 1), c(14, 2), c(14, 3), c(14, 4), c(17, 6), c(17, 7),
  c(18, 1), c(18, 2), c(20, 1), c(20, 2), c(22, 1), c(22, 2), c(26, 24),
  c(26, 25), c(28, 3), c(28, 24), c(28, 25), c(29, 3), c(30, 24),
  c(30, 27), c(31, 2), c(31, 9), c(32, 1), c(32, 25), c(32, 26),
  c(32, 29), c(33, 3), c(33, 9), c(33, 15), c(33, 16), c(33, 19),
  c(33, 21), c(33, 23), c(33, 24), c(33, 30), c(33, 31), c(33, 32),
  c(34, 9), c(34, 10), c(34, 14), c(34, 15), c(34, 16), c(34, 19),
  c(34, 20), c(34, 21), c(34, 23), c(34, 24), c(34, 27), c(34, 28),
  c(34, 29), c(34, 30), c(34, 31), c(34, 32), c(34, 33)
)

# Populate the matrix with 1s based on connections
for (connection in connections) {
  # (I'm transposing them so that I am filling out the upper triangular portion
  #  rather than the lower triangular portion, but this doesn't matter too much.)
  zachary[connection[1], connection[2]] <- 1
}

# Set the upper diagonal of the matrix to also be the same
zachary <- zachary + t(zachary)

# "True community membership" as determined by "club after fission"
# from Table 1 in Wayne Zachary's paper.
# Here, a "1" denotes membership in Mr. Hi's faction and a "2" is membership
# in Officers' faction.
true_communities <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2,
                      1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

# ============================================
# == Check through multiple values of gamma ==
# ============================================

gamma_check <- seq(0.00, 0.50, length.out = 50)
num_sim_per_gamma <- 200
rand_results_fission_true <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
rand_results_fission_full <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_midpoint <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_lower <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_upper <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]

  for (rep in 1:num_sim_per_gamma) {

    # Do clustering on the full version of the matrix
    z_hat_full_initial <- nett::spec_clust(zachary, K = K)
    z_hat_full <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat_full[, i] <- 1 * (z_hat_full_initial == i)
    }
    n_hat_full <- apply(z_hat_full, 2, sum)

    # -------------------------------------------------
    # -- Fission the matrix with this value of gamma --
    # -------------------------------------------------

    # Although we create an entire 34x34 matrix, only the upper triangular
    # portion ever gets used.
    Zfission_temp <- matrix(rbinom(n = n^2, size = 1, prob = gamma), nrow = n)
    Zfission <- matrix(0, nrow = n, ncol = n)
    Zfission[upper.tri(Zfission)] <- Zfission_temp[upper.tri(Zfission_temp)]
    Zfission <- Zfission + t(Zfission)

    # Fission
    zachary_tr <- zachary * (1 - Zfission) + (1 - zachary) * Zfission

    # Cluster using the training set
    z_hat_initial <- nett::spec_clust(zachary_tr, K = K)
    z_hat <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat[, i] <- 1 * (z_hat_initial == i)
    }
    n_hat <- apply(z_hat, 2, sum)

    # -----------------------------------
    # -- Check agreement of clustering --
    # -----------------------------------

    rand_results_fission_full[gamma_index, rep] <- mclust::adjustedRandIndex(z_hat_full_initial, z_hat_initial)
    rand_results_fission_true[gamma_index, rep] <- mclust::adjustedRandIndex(true_communities, z_hat_initial)

    # ------------------
    # -- Do inference --
    # ------------------

    # Calculate our estimator and its target
    B <- matrix(0, nrow = K, ncol = K)
    Delta <- matrix(0, nrow = K, ncol = K)

    for (k in 1:K) {
      for (l in 1:K) {
        k_nodes <- which(z_hat_initial == k)
        l_nodes <- which(z_hat_initial == l)

        # Create the index set of all edges which point from k to l,
        # and all edges that point from l to k.
        Ikl <- rbind(as.matrix(expand.grid(k_nodes, l_nodes)),
                     as.matrix(expand.grid(l_nodes, k_nodes)))
        # Cut down this index set to just the edges (i, j) where i < j
        # which avoids the redundancy due to the symmetry in the network.
        Ikl_prime <- Ikl[Ikl[, 1] < Ikl[, 2], ]

        # Estimator
        for (ij_index in 1:nrow(Ikl_prime)) {
          ij <- Ikl_prime[ij_index, ]
          B[k, l] <- B[k, l] + zachary[ij[1], ij[2]]
        }
        B[k, l] <- B[k, l] / nrow(Ikl_prime)

        # Variance
        Delta[k, l] <- (B[k, l] * (1 - B[k, l])) / nrow(Ikl_prime)
      }
    }

    u <- c(1, -2, 0, 1)
    u <- u / sqrt(sum(u^2))

    # Build a confidence interval
    alpha <- 0.10
    estimate <- t(u) %*% as.vector(B)
    estimate_var <- t(u) %*% diag(as.vector(Delta)) %*% u
    margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)

    # Record the confidence interval
    ci_midpoint[gamma_index, rep] <- estimate
    ci_lower[gamma_index, rep] <- estimate - margin_of_error
    ci_upper[gamma_index, rep] <- estimate + margin_of_error
  }
}

# ====================
# == Create figures ==
# ====================

# ---------------------------------------
# -- RAND index as a function of gamma --
# ---------------------------------------

avg_rand_results_fission_full <- apply(rand_results_fission_full, 1, mean)
avg_rand_results_fission_true <- apply(rand_results_fission_true, 1, mean)

plot_df <- data.frame(gamma = gamma_check,
                      avg_rand_fission_full = avg_rand_results_fission_full,
                      avg_rand_fission_true = avg_rand_results_fission_true)

zachary_randindex_gamma <- ggplot(plot_df) +
  geom_line(aes(x = gamma, y = avg_rand_fission_full, color = 'against_full'), alpha = 0.7, size = 1.1) +
  geom_line(aes(x = gamma, y = avg_rand_fission_true, color = 'against_true'), alpha = 0.7, size = 1.1) +
  xlab(TeX('$\\gamma$')) + ylab('Adjusted RAND Index') +
  labs(color = 'Comparison') +
  scale_color_manual(
    values = c('against_full' = 'darkslategray4', against_true = 'firebrick3'),
    labels = c(TeX('$\\hat{Z}^{(tr),\\gamma}$ and $\\hat{Z}^{(full)}$'), TeX('$\\hat{Z}^{(tr),\\gamma}$ and $Z^{(true)}$'))
  ) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
zachary_randindex_gamma
ggsave('figures/zachary_randindex_gamma.pdf', plot = zachary_randindex_gamma,
       device = 'pdf', width = 6, height = 4.5)

# -------------------------------------------------
# -- Confidence intervals as a function of gamma --
# -------------------------------------------------

avg_ci_midpoint <- apply(ci_midpoint, 1, mean)
avg_ci_lower <- apply(ci_lower, 1, mean)
avg_ci_upper <- apply(ci_upper, 1, mean)

plot_df <- data.frame(gamma = gamma_check,
                      avg_ci_midpoint = avg_ci_midpoint,
                      avg_ci_lower = avg_ci_lower,
                      avg_ci_upper = avg_ci_upper)

zachary_CI_gamma <- ggplot(plot_df) +
  geom_line(aes(x = gamma, y = avg_ci_lower, color = 'bounds'), size = 1.1) +
  geom_line(aes(x = gamma, y = avg_ci_midpoint, color = 'midpoint'), size = 1.1) +
  geom_line(aes(x = gamma, y = avg_ci_upper, color = 'bounds'), size = 1.1) +
  geom_line(aes(x = gamma, y = 0), linetype = 'dashed', alpha = 0.6, size = 0.8) +
  xlab(TeX('$\\gamma$')) + ylab(TeX('$\\theta(A^{(tr)}_\\gamma)$')) +
  labs(color = 'Confidence interval') +
  scale_color_manual(
    values = c('bounds' = 'darkslategray4', 'midpoint' = 'firebrick3'),
    labels = c('Lower/upper bounds', 'Estimate')
  ) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
zachary_CI_gamma
ggsave('figures/zachary_CI_gamma.pdf', plot = zachary_CI_gamma,
       device = 'pdf', width = 6, height = 4.5)
