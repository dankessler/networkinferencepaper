library('networkinference')
library('latex2exp')
library('tidyverse')
library('igraph')

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
zachary_upper_diagonal <- zachary
zachary <- zachary + t(zachary)

# "True community membership" as determined by "club after fission"
# from Table 1 in Wayne Zachary's paper.
# Here, a "1" denotes membership in Mr. Hi's faction and a "2" is membership
# in Officers' faction.
true_communities <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2,
                      1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

# =====================================================
# == Create a cartoon for Zachary's karate club data ==
# =====================================================

# Visualize adjaceny matrix
zachary_visual_df <- reshape2::melt(zachary)
colnames(zachary_visual_df) <- c('Row', 'Column', 'Value')

# With numbers for all the cells
zachary_visual_2 <- ggplot(zachary_visual_df) +
  geom_tile(aes(x = Column, y = Row, fill = factor(Value)), color = 'gray80') +
  scale_fill_manual(values = c('1' = 'black', '0' = 'white')) +
  scale_y_reverse(breaks = 1:n, labels = 1:n, expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:n, labels = 1:n, expand = c(0, 0), position = 'top') +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(size = 8),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = NULL)

# Without numbers
zachary_visual_2 <- ggplot(zachary_visual_df) +
  geom_tile(aes(x = Column, y = Row, fill = factor(Value)), color = 'gray80') +
  scale_fill_manual(values = c('1' = 'black', '0' = 'white')) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), position = 'top') +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = NULL)

zachary_visual_2
ggsave('figures/zachary/zachary_cartoon_2.pdf', zachary_visual_2, device = 'pdf',
       width = 7, height = 7)

# More direct way of showing network
set.seed(6)
zach_viz <- graph_from_adjacency_matrix(zachary, mode = c('undirected'),
                                        weighted = FALSE, diag = FALSE)
zach_viz_plot <- plot(zach_viz, vertex.size = 13, vertex.color = 'seashell1',
                      vertex.label.font = 2, vertex.label.color = 'black',
                      vertex.label.family = 'Helvetica',
                      vertex.label.cex = 1.4, edge.width=1.2, edge.color='black', margin = 0)


# =======================
# == Helpful functions ==
# =======================

logit <- function(x) {
  return(log(x / (1 - x)))
}
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}
h0 <- function(x, gamma) {
  return(expit(logit(x) + log(gamma / (1 - gamma))))
}
h0_alt <- function(x, gamma) {
  return(x / (x + (1-x)*((1-gamma) / gamma)))
}
h1 <- function(x, gamma) {
  return(expit(logit(x) + log((1 - gamma) / gamma)))
}
h0_inv <- function(x, gamma) {
  return(h1(x, gamma))
}
h1_inv <- function(x, gamma) {
  return(h0(x, gamma))
}
h1_inv_deriv <- function(x, gamma) {
  c0 <- log(gamma / (1 - gamma))
  return((expit(logit(x) + c0) / (1 + exp(logit(x) + c0))) / (x*(1-x)))
}
h0_inv_deriv <- function(x, gamma) {
  c1 <- log((1 - gamma) / gamma)
  return((expit(logit(x) + c1) / (1 + exp(logit(x) + c1))) / (x*(1-x)))
}

# ============================================
# == Check through multiple values of gamma ==
# ============================================

gamma_check <- seq(0.0001, 0.50, length.out = 20)
gamma_check_failed <- rep(FALSE, length.out = length(gamma_check))
num_sim_per_gamma <- 500

rand_results_fission_true <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
rand_results_fission_full <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

ci_midpoint <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_lower <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_upper <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]
  c0 <- log(gamma / (1-gamma))
  c1 <- log((1-gamma) / gamma)

  for (rep in 1:num_sim_per_gamma) {

    # Sometimes, iterations of fission result in cases where we have all zeros
    # or ones within a set of nodes where A_tr = 0 or A_tr = 1, which causes
    # division by zero issues. For purposes of "averaging over gamma," I will
    # exclude these cases since it is unclear how to handle them.
    need_good_matrix <- TRUE
    times_good_matrix_attempted <- 0
    limits_on_attempts <- 1000

    while(need_good_matrix) {
      times_good_matrix_attempted <- times_good_matrix_attempted + 1

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
      B0 <- matrix(0, nrow = K, ncol = K)
      B1 <- matrix(0, nrow = K, ncol = K)
      Delta0 <- matrix(0, nrow = K, ncol = K)
      Delta1 <- matrix(0, nrow = K, ncol = K)
      Delta <- matrix(0, nrow = K, ncol = K)
      Phi <- matrix(0, nrow = K, ncol = K)

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
          num_Ikl_pr <- nrow(Ikl_prime)
          num_0s <- 0
          num_1s <- 0

          # Estimator
          for (ij_index in 1:nrow(Ikl_prime)) {
            ij <- Ikl_prime[ij_index, ]

            if (zachary_tr[ij[1], ij[2]] == 0) {
              B0[k, l] <- B0[k, l] + zachary[ij[1], ij[2]]
              num_0s <- num_0s + 1
            } else {
              B1[k, l] <- B1[k, l] + zachary[ij[1], ij[2]]
              num_1s <- num_1s + 1
            }
          }
          B0[k, l] <- B0[k, l] / num_0s
          B1[k, l] <- B1[k, l] / num_1s

          # Only factor in Ikl_s if num_s > 0
          Phi_kl_0 <- 0
          Phi_kl_1 <- 0
          Delta_kl_0 <- 0
          Delta_kl_1 <- 0
          if (num_0s > 0) {
            # Mean estimate
            Phi_kl_0 <- (num_0s / num_Ikl_pr) * (B0[k, l] / (B0[k, l] + ((1 - B0[k, l]) * gamma / (1-gamma))))

            # Variance estimate
            B0_kl_adjust <- B0[k, l]
            if (B0[k, l] == 0) {
              B0_kl_adjust <- 1 / (2*num_0s)
            }
            if (B0[k, l] == 1) {
              B0_kl_adjust <- (num_0s - 0.5) / num_0s
            }
            Delta_kl_0 <- (B0_kl_adjust * (1 - B0_kl_adjust) * exp(2*c0)) /
              (num_0s * ((1 - B0[k, l])*exp(c0) + B0[k, l])^4)
          }
          if (num_1s > 0) {
            # Mean estimate
            Phi_kl_1 <- (num_1s / num_Ikl_pr) * (B1[k, l] / (B1[k, l] + ((1 - B1[k, l]) * (1-gamma) / gamma)))

            # Variance estimate
            B1_kl_adjust <- B1[k, l]
            if (B1[k, l] == 0) {
              B1_kl_adjust <- 1 / (2*num_1s)
            }
            if (B1[k, l] == 1) {
              B1_kl_adjust <- (num_1s - 0.5) / num_1s
            }
            Delta_kl_1 <- (B1_kl_adjust * (1 - B1_kl_adjust) * exp(2*c1)) /
              (num_1s * ((1 - B1[k, l])*exp(c1) + B1[k, l])^4)
          }

          # Filter the Delta_kl_0 and Delta_kl_1 individually to be bounded
          # by 0.25
          Delta_kl_0 <- pmin(Delta_kl_0, 0.25)
          Delta_kl_1 <- pmin(Delta_kl_1, 0.25)

          # Combine the 0s and 1s to yield the final mean estimate
          # as well as the final variance estimate
          Phi[k, l] <- Phi_kl_0 + Phi_kl_1
          Delta[k, l] <- (num_0s / num_Ikl_pr)^2 * Delta_kl_0 +
            (num_1s / num_Ikl_pr)^2 * Delta_kl_1


          # For the variance calculation only, if B0 or B1 is exactly 0 or 1, then
          # shrink it a small amount toward 1/2.
          # B0_kl_var <- B0[k, l]
          # B1_kl_var <- B1[k, l]
          # if (B0[k, l] == 0) {
          #   B0_kl_var <- 1 / num_0s
          # }
          # if (B1[k, l] == 1) {
          #   B1_kl_var <- (num_1s - 1) / num_1s
          # }
          #
          # # Variance
          # G1_G2_0 <- (B0_kl_var * (1 - B0_kl_var) * exp(2*c0)) /
          #   (num_0s * ((1 - B0_kl_var)*exp(c0) + B0_kl_var)^4)
          # G1_G2_1 <- (B1_kl_var * (1 - B1_kl_var) * exp(2*c1)) /
          #   (num_1s * ((1 - B1_kl_var)*exp(c1) + B1_kl_var)^4)
          #
          # G1_G2_0[is.nan(G1_G2_0)] <- Inf
          # G1_G2_1[is.nan(G1_G2_1)] <- Inf

          # Should we make really wide confidence intervals when G1_G2_0 is 0?
          # G1_G2_0[G1_G2_0 == 0] <- Inf
          # G1_G2_1[G1_G2_1 == 0] <- Inf

          # G1_G2_0 <- pmin(G1_G2_0, 0.25)
          # G1_G2_1 <- pmin(G1_G2_1, 0.25)


          # Delta0 <- ((B0[k, l] * (1 - B0[k, l])) / num_0s) * (h0_inv_deriv(B0[k, l], gamma))^2
          # Delta1 <- ((B1[k, l] * (1 - B1[k, l])) / num_1s) * (h1_inv_deriv(B1[k, l], gamma))^2
          # Delta[k, l] <- (num_0s / num_Ikl_pr)^2 * G1_G2_0 + (num_1s / num_Ikl_pr)^2 * G1_G2_1
        }
      }

      # Check if we passed the variance check
      if (sum(is.nan(Delta)) == 0) {
        need_good_matrix <- FALSE

        u <- c(1, -2, 0, 1)
        u <- u / sqrt(sum(u^2))

        # Build a confidence interval
        alpha <- 0.90
        estimate <- t(u) %*% as.vector(Phi)
        estimate_var <- t(u) %*% diag(as.vector(Delta)) %*% u
        margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)

        # Record the confidence interval
        ci_midpoint[gamma_index, rep] <- estimate
        ci_lower[gamma_index, rep] <- estimate - margin_of_error
        ci_upper[gamma_index, rep] <- estimate + margin_of_error
      }
      if (times_good_matrix_attempted > limits_on_attempts) {
        need_good_matrix <- FALSE
        gamma_check_failed[gamma_index] <- TRUE
      }
    }
  }
}

# ====================================================
# == Also do estimation with true known communities ==
# ====================================================

u <- c(1, -2, 0, 1)
u <- u / sqrt(sum(u^2))

z_true <- matrix(rep(NA, n*K), nrow = n)
for (i in 1:K) {
  z_true[, i] <- 1 * (true_communities == i)
}
n_true <- apply(z_true, 2, sum)
NN_true_inv <- diag(1 / diag(t(z_true) %*% z_true))

zachary_matrix_est <- matrix(0, nrow = K, ncol = K)
zachary_matrix_var_est <- matrix(0, nrow = K, ncol = K)

for (k in 1:K) {
  for (l in 1:K) {
    k_nodes <- which(true_communities == k)
    l_nodes <- which(true_communities == l)

    # Create the index set of all edges which point from k to l,
    # and all edges that point from l to k.
    Ikl <- rbind(as.matrix(expand.grid(k_nodes, l_nodes)),
                 as.matrix(expand.grid(l_nodes, k_nodes)))
    # Cut down this index set to just the edges (i, j) where i < j
    # which avoids the redundancy due to the symmetry in the network.
    Ikl_prime <- Ikl[Ikl[, 1] < Ikl[, 2], ]
    num_Ikl_pr <- nrow(Ikl_prime)

    # Estimator
    for (ij_index in 1:nrow(Ikl_prime)) {
      ij <- Ikl_prime[ij_index, ]
      zachary_matrix_est[k, l] <- zachary_matrix_est[k, l] + zachary[ij[1], ij[2]]
    }
    zachary_matrix_est[k, l] <- zachary_matrix_est[k, l] / num_Ikl_pr
    zachary_matrix_var_est[k, l] <- zachary_matrix_est[k, l] * (1 - zachary_matrix_est[k, l]) / num_Ikl_pr
  }
}

# Build a confidence interval
alpha <- 0.90
estimate_true <- t(u) %*% as.vector(zachary_matrix_est)
estimate_true_var <- t(u) %*% diag(as.vector(zachary_matrix_var_est)) %*% u
margin_of_error_true <- qnorm(1 - alpha / 2) * sqrt(estimate_true_var)

estimate_true_lb <- estimate_true - margin_of_error_true
estimate_true_ub <- estimate_true + margin_of_error_true

# Z(A) vs Ztrue
mclust::adjustedRandIndex(z_hat_full_initial, true_communities)

# For 1/30/25 suggestion I look at theta(A) using Zhat
theta_Zhat_midpoint <- matrix(0, nrow = K, ncol = K)
for (k in 1:K) {
  for (l in 1:K) {
    k_nodes <- which(z_hat_full_initial == k)
    l_nodes <- which(z_hat_full_initial == l)

    # Create the index set of all edges which point from k to l,
    # and all edges that point from l to k.
    Ikl <- rbind(as.matrix(expand.grid(k_nodes, l_nodes)),
                 as.matrix(expand.grid(l_nodes, k_nodes)))
    # Cut down this index set to just the edges (i, j) where i < j
    # which avoids the redundancy due to the symmetry in the network.
    Ikl_prime <- Ikl[Ikl[, 1] < Ikl[, 2], ]
    num_Ikl_pr <- nrow(Ikl_prime)

    # Estimator
    for (ij_index in 1:nrow(Ikl_prime)) {
      ij <- Ikl_prime[ij_index, ]
      theta_Zhat_midpoint[k, l] <- theta_Zhat_midpoint[k, l] + zachary[ij[1], ij[2]]
    }
    theta_Zhat_midpoint[k, l] <- theta_Zhat_midpoint[k, l] / num_Ikl_pr
  }
}
theta_Zhat_temp <- t(u) %*% as.vector(theta_Zhat_midpoint)

# ====================
# == Create figures ==
# ====================

axis.text.x.size <- 15
axis.text.y.size <- 15
axis.title.x.size <- 21
axis.title.y.size <- 14
legend.title.size <- 18
legend.text.size <- 15
legend.text.size.small <- 14
thick_linewidth = 1.0
thin_linewidth = 0.9
legend.linewidth.thick <- 0.9
legend.linewidth.thin <- 0.7

# ---------------------------------------
# -- RAND index as a function of gamma --
# ---------------------------------------

avg_rand_results_fission_full <- apply(rand_results_fission_full, 1, mean)
avg_rand_results_fission_true <- apply(rand_results_fission_true, 1, mean)

plot_df <- data.frame(gamma = gamma_check,
                      avg_rand_fission_full = avg_rand_results_fission_full,
                      avg_rand_fission_true = avg_rand_results_fission_true)

zachary_randindex_gamma <- ggplot(plot_df) +
  geom_line(aes(x = gamma, y = avg_rand_fission_full, color = 'against_full'), alpha = 0.9, linewidth = thick_linewidth) +
  geom_line(aes(x = gamma, y = avg_rand_fission_true, color = 'against_true'), alpha = 0.9, linewidth = thick_linewidth) +
  xlab(TeX('$\\gamma$')) + ylab('Adjusted RAND Index') +
  labs(color = 'Legend') +
  scale_color_manual(
    values = c('against_full' = 'darkslategray4', against_true = 'firebrick3'),
    labels = c(TeX('$\\hat{Z}(A^{(tr)}_\\gamma)$ and $\\hat{Z}(A)$'), TeX('$\\hat{Z}(A^{(tr)}_\\gamma)$ and $Z^{(true)}$'))
  ) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size))
zachary_randindex_gamma
ggsave('figures/zachary_randindex_gamma.pdf', plot = zachary_randindex_gamma,
       device = 'pdf', width = 14, height = 8, units = 'cm')

# -------------------------------------------------
# -- Confidence intervals as a function of gamma --
# -------------------------------------------------

axis.text.x.size <- 15
axis.text.y.size <- 15
axis.title.x.size <- 21
axis.title.y.size <- 14
legend.title.size <- 18
legend.text.size <- 15
legend.text.size.small <- 14
thick_linewidth = 1.0
thin_linewidth = 0.7
legend.linewidth.thick <- 0.9
legend.linewidth.thin <- 0.7

avg_ci_midpoint <- apply(ci_midpoint, 1, mean, na.rm = TRUE)
avg_ci_lower <- apply(ci_lower, 1, mean, na.rm = TRUE)
avg_ci_upper <- apply(ci_upper, 1, mean, na.rm = TRUE)

plot_df_midpoint <- data.frame(gamma = gamma_check, avg_xi_hat = avg_ci_midpoint)
plot_df_ci_lower <- data.frame(gamma = gamma_check, ci_bounds_xi = avg_ci_lower)
plot_df_ci_upper <- data.frame(gamma = gamma_check, ci_bounds_xi = avg_ci_upper)

plot_df_errorbars <- data.frame(gamma = gamma_check,
                                avg_xi_hat = avg_ci_midpoint,
                                lowpoint = avg_ci_lower,
                                highpoint = avg_ci_upper)
plot_df_errorbars_true <- data.frame(gamma = gamma_check,
                                     est_true = rep(estimate_true, length(gamma_check)),
                                     est_true_lb = rep(estimate_true_lb, length(gamma_check)),
                                     est_true_ub = rep(estimate_true_ub, length(gamma_check)))

# Workaround for getting long math label to get on two lines
# theta_true_bounds_top_exp <- TeX('90% CI for $\\theta(A)$')
# theta_true_bounds_bottom_exp <- TeX('using $Z^{(true)}$')
# theta_true_bounds_exp <- expression(atop())

legend_values <- c('avg_xi_hat' = 'palevioletred3',
                   'est_true' = 'cadetblue',
                   'bounds' = 'black')
legend_values_linetype <- c('avg_xi_hat' = 'solid',
                            'est_true' = 'solid')
legend_labels <- c('avg_xi_hat' = TeX('$\\hat{\\xi}(A_{\\gamma}^{(tr)})$'),
                   'est_true' = TeX('$\\hat{\\theta}(A)$ using $Z^{(true)}$'))

zachary_CI_gamma <- ggplot() +
  geom_line(data = plot_df_midpoint, aes(x = gamma, y = 0), linetype = 'dashed', alpha = 0.6, linewidth = thin_linewidth) +
  geom_line(data = plot_df_midpoint, aes(x = gamma, y = avg_xi_hat, color = 'avg_xi_hat', linetype = 'avg_xi_hat'), linewidth = thick_linewidth) +
  geom_line(data = plot_df_midpoint, aes(x = gamma, y = estimate_true, color = 'est_true', linetype = 'est_true'), linewidth = thick_linewidth) +
  geom_linerange(data = plot_df_errorbars_true, aes(x = gamma, y = est_true, ymin = est_true_lb, ymax = est_true_ub, linetype = 'est_true'), color = 'black', linewidth = thin_linewidth, show.legend = FALSE) +
  geom_linerange(data = plot_df_errorbars, aes(x = gamma, y = avg_xi_hat, ymin = lowpoint, ymax = highpoint, linetype = 'est_true'), color = 'black', linewidth = thin_linewidth, show.legend = FALSE) +
  xlab(TeX('$\\gamma$')) + ylab('') +
  labs(color = 'Legend') +
  scale_color_manual(
    values = legend_values,
    labels = legend_labels
  ) +
  scale_linetype_manual(
    values = legend_values_linetype,
    labels = legend_labels
  ) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size)) +
  guides(linetype = 'none')
zachary_CI_gamma
ggsave('figures/zachary_CI_gamma.pdf', plot = zachary_CI_gamma,
       device = 'pdf', width = 14, height = 8, units = 'cm')

# Old one with dashed lines for the confidence intervals

# # -------------------------------------------------
# # -- Confidence intervals as a function of gamma --
# # -------------------------------------------------
#
# axis.text.x.size <- 15
# axis.text.y.size <- 15
# axis.title.x.size <- 21
# axis.title.y.size <- 14
# legend.title.size <- 18
# legend.text.size <- 15
# legend.text.size.small <- 14
# thick_linewidth = 1.0
# thin_linewidth = 0.5
# legend.linewidth.thick <- 0.9
# legend.linewidth.thin <- 0.7
#
# avg_ci_midpoint <- apply(ci_midpoint, 1, mean, na.rm = TRUE)
# avg_ci_lower <- apply(ci_lower, 1, mean, na.rm = TRUE)
# avg_ci_upper <- apply(ci_upper, 1, mean, na.rm = TRUE)
#
# plot_df_midpoint <- data.frame(gamma = gamma_check, avg_xi_hat = avg_ci_midpoint)
# plot_df_ci_lower <- data.frame(gamma = gamma_check, ci_bounds_xi = avg_ci_lower)
# plot_df_ci_upper <- data.frame(gamma = gamma_check, ci_bounds_xi = avg_ci_upper)
#
# plot_df_errorbars <- data.frame(gamma = gamma_check,
#                                midpoint = avg_ci_midpoint,
#                                lowpoint = avg_ci_lower,
#                                highpoint = avg_ci_upper)
#
# # Workaround for getting long math label to get on two lines
# # theta_true_bounds_top_exp <- TeX('90% CI for $\\theta(A)$')
# # theta_true_bounds_bottom_exp <- TeX('using $Z^{(true)}$')
# # theta_true_bounds_exp <- expression(atop())
#
# legend_values <- c('avg_xi_hat' = 'palevioletred3',
#                    'ci_bounds_xi' = 'palevioletred3',
#                    'est_true' = 'cadetblue',
#                    'est_true_bounds' = 'cadetblue')
# legend_values_linetype <- c('avg_xi_hat' = 'solid',
#                             'ci_bounds_xi' = 'dashed',
#                             'est_true' = 'solid',
#                             'est_true_bounds' = 'dashed')
# legend_labels <- c('avg_xi_hat' = TeX('$\\hat{\\xi}(A_{\\gamma}^{(tr)})$'),
#                    'ci_bounds_xi' = TeX('90% CI for $\\xi(A_{\\gamma}^{(tr)})$'),
#                    'est_true' = TeX('$\\hat{\\theta}(A)$ using $Z^{(true)}$'),
#                    'est_true_bounds' = TeX('90% CI for $\\theta(A)$'))
#                    # 'est_true_bounds' = TeX('90% CI for $\\theta(A)$ using $Z^{(true)}$'))
#
# zachary_CI_gamma <- ggplot() +
#   geom_line(data = plot_df_midpoint, aes(x = gamma, y = 0), linetype = 'dashed', alpha = 0.6, linewidth = thin_linewidth) +
#   geom_line(data = plot_df_midpoint, aes(x = gamma, y = avg_xi_hat, color = 'avg_xi_hat', linetype = 'avg_xi_hat'), linewidth = thick_linewidth) +
#   geom_line(data = plot_df_ci_upper, aes(x = gamma, y = ci_bounds_xi, color = 'ci_bounds_xi', linetype = 'ci_bounds_xi'), linewidth = thin_linewidth) +
#   geom_line(data = plot_df_ci_lower, aes(x = gamma, y = ci_bounds_xi, color = 'ci_bounds_xi', linetype = 'ci_bounds_xi'), linewidth = thin_linewidth) +
#   geom_line(data = plot_df_midpoint, aes(x = gamma, y = estimate_true, color = 'est_true', linetype = 'est_true'), linewidth = thick_linewidth) +
#   geom_line(data = plot_df_midpoint, aes(x = gamma, y = estimate_true_lb, color = 'est_true_bounds', linetype = 'est_true_bounds'), linewidth = thin_linewidth) +
#   geom_line(data = plot_df_midpoint, aes(x = gamma, y = estimate_true_ub, color = 'est_true_bounds', linetype = 'est_true_bounds'), linewidth = thin_linewidth) +
#   geom_linerange(data = plot_df_errorbars, aes(x = gamma, y = midpoint, ymin = lowpoint, ymax = highpoint)) +
#   xlab(TeX('$\\gamma$')) + ylab('') +
#   labs(color = 'Legend') +
#   scale_color_manual(
#     values = legend_values,
#     labels = legend_labels
#   ) +
#   scale_linetype_manual(
#     values = legend_values_linetype,
#     labels = legend_labels
#   ) +
#   theme(aspect.ratio = 1,
#         axis.text.x = element_text(size = axis.text.x.size),
#         axis.text.y = element_text(size = axis.text.y.size),
#         axis.title.x = element_text(size = axis.title.x.size),
#         axis.title.y = element_text(size = axis.title.y.size),
#         legend.title = element_text(size = legend.title.size),
#         legend.text = element_text(size = legend.text.size)) +
#   guides(color = guide_legend(override.aes = list(
#     linetype = c('solid', 'dashed'),
#     linewidth = c(legend.linewidth.thick, legend.linewidth.thin)
#   )),
#   linetype = 'none')
# zachary_CI_gamma
# ggsave('figures/zachary_CI_gamma.pdf', plot = zachary_CI_gamma,
#        device = 'pdf', width = 14, height = 8, units = 'cm')
#
#
# ci_widths <- avg_ci_upper - avg_ci_lower
# plot_df <- data.frame(gamma = gamma_check, ci_widths = ci_widths)
# ggplot(plot_df) +
#   geom_line(aes(x = gamma, y = ci_widths))
