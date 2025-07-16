library('networkinference')
library('latex2exp')
library('tidyverse')
library('igraph')

ggplot2::theme_set(theme_minimal())

# == Dolphins? ==
library("manynet")
data(ison_dolphins)

# Nodes
dolphins_igraph <- as_igraph(ison_dolphins)
dolphins <- igraph::as_adjacency_matrix(dolphins_igraph, sparse = FALSE)

# =====================================================
# == Manually load the Zachary's karate club dataset ==
# =====================================================

n <- 62
# n <- 34
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
  zachary[connection[2], connection[1]] <- 1
}

# Make the matrix symmetric
# (not necessary for networkinference package, but necessary for spectral clustering.)
zachary[lower.tri(zachary)] <- t(zachary)[lower.tri(zachary)]

# "True community membership" as determined by "club after fission"
# from Table 1 in Wayne Zachary's paper.
# Here, a "1" denotes membership in Mr. Hi's faction and a "2" is membership
# in Officers' faction.
# true_communities <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2,
#                       1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

# ========================================================== #
# == Conduct inference for selected target, varying gamma == #
# ========================================================== #

gamma_check <- seq(0.01, 0.50, length.out = 20)
num_sim_per_gamma <- 200

rand_results_fission_true <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
rand_results_fission_full <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

ci_midpoint <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_lower <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
ci_upper <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

for (gamma_index in 1:length(gamma_check)) {
  gamma <- gamma_check[gamma_index]
  c0 <- log(gamma / (1-gamma))
  c1 <- log((1-gamma) / gamma)
  cat(paste0("Checking gamma = ", gamma, "\n"))

  for (rep in 1:num_sim_per_gamma) {

    # Linear combination vector
    u <- c(1, 0, -2, 1)
    u <- u / sqrt(sum(u^2))

    # Do clustering on the full version of the matrix
    z_hat_full_initial <- nett::spec_clust(dolphins, K = K)
    z_hat_full <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat_full[, i] <- 1 * (z_hat_full_initial == i)
    }
    n_hat_full <- apply(z_hat_full, 2, sum)

    # Split matrix
    dolphins_splitting <- networkinference::split_matrix(A = dolphins,
                                                        distribution = "bernoulli",
                                                        gamma = gamma,
                                                        allow_self_loops = FALSE,
                                                        is_directed = FALSE)
    dolphins_tr <- dolphins_splitting$Atr
    dolphins_te <- dolphins_splitting$Ate

    # Make the train and test sets symmetric (makes spectral clustering work out better
    # when implemented in the nett package it seems)
    dolphins_tr[lower.tri(dolphins_tr)] <- t(dolphins_tr)[lower.tri(dolphins_tr)]
    dolphins_te[lower.tri(dolphins_te)] <- t(dolphins_te)[lower.tri(dolphins_te)]

    # Cluster using the training set
    z_hat_initial <- nett::spec_clust(dolphins_tr, K = K)
    z_hat <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat[, i] <- 1 * (z_hat_initial == i)
    }
    n_hat <- apply(z_hat, 2, sum)

    # Skip this iteration if there are no nodes assigned to one of the two
    # communities
    if (any(n_hat == 1) | any(n_hat == 0)) {
      cat(paste0("\tSkipping iteration ", rep, ". One of the estimated communities contains only a single node.\n"))
      rep <- rep - 1
      next
    }

    # Check agreement of clustering
    rand_results_fission_full[gamma_index, rep] <- mclust::adjustedRandIndex(z_hat_full_initial, z_hat_initial)
    # rand_results_fission_true[gamma_index, rep] <- mclust::adjustedRandIndex(true_communities, z_hat_initial)

    # Conduct inference for selected target
    dolphins_inference <- networkinference::infer_network(Ate = dolphins_te,
                                                         u = u,
                                                         communities = z_hat_initial,
                                                         distribution = "bernoulli",
                                                         K = 2, gamma = gamma,
                                                         Atr = dolphins_tr,
                                                         allow_self_loops = FALSE,
                                                         is_directed = FALSE)

    # Build a confidence interval
    alpha <- 0.10
    estimate <- dolphins_inference$estimate
    estimate_var <- dolphins_inference$estimate_variance
    margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)

    # Record the confidence interval
    ci_midpoint[gamma_index, rep] <- estimate
    ci_lower[gamma_index, rep] <- estimate - margin_of_error
    ci_upper[gamma_index, rep] <- estimate + margin_of_error
  }
}

# ================================================ #
# == Also do estimation with "true" communities == #
# ================================================ #

# Note that this is manually implemented because my package does not
# accommodate this type of procedure (where true communities are known)

u <- c(1, 0, -2, 1)
u <- u / sqrt(sum(u^2))

z_true <- matrix(rep(NA, n*K), nrow = n)
for (i in 1:K) {
  z_true[, i] <- 1 * (true_communities == i)
}
n_true <- apply(z_true, 2, sum)
NN_true_inv <- diag(1 / diag(t(z_true) %*% z_true))

dolphins_matrix_est <- matrix(0, nrow = K, ncol = K)
dolphins_matrix_var_est <- matrix(0, nrow = K, ncol = K)

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
      dolphins_matrix_est[k, l] <- dolphins_matrix_est[k, l] + dolphins[ij[1], ij[2]]
    }
    dolphins_matrix_est[k, l] <- dolphins_matrix_est[k, l] / num_Ikl_pr
    dolphins_matrix_var_est[k, l] <- dolphins_matrix_est[k, l] * (1 - dolphins_matrix_est[k, l]) / num_Ikl_pr
  }
}

# Build a confidence interval
alpha <- 0.10
estimate_true <- t(u) %*% as.vector(dolphins_matrix_est)
estimate_true_var <- t(u) %*% diag(as.vector(dolphins_matrix_var_est)) %*% u
margin_of_error_true <- qnorm(1 - alpha / 2) * sqrt(estimate_true_var)

estimate_true_lb <- estimate_true - margin_of_error_true
estimate_true_ub <- estimate_true + margin_of_error_true

# ==================== #
# == Create figures == #
# ==================== #

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

# ---------------------------------------------------- #
# -- Left panel - RAND index as a function of gamma -- #
# ---------------------------------------------------- #

avg_rand_results_fission_full <- apply(rand_results_fission_full, 1, mean)
avg_rand_results_fission_true <- apply(rand_results_fission_true, 1, mean)

plot_df <- data.frame(gamma = gamma_check,
                      avg_rand_fission_full = avg_rand_results_fission_full,
                      avg_rand_fission_true = avg_rand_results_fission_true)

dolphins_randindex_gamma <- ggplot(plot_df) +
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
dolphins_randindex_gamma
ggsave('figures/dolphins_randindex_gamma.pdf', plot = dolphins_randindex_gamma,
       device = 'pdf', width = 14, height = 8, units = 'cm')

# --------------------------------------------------------------- #
# -- Right panel - Confidence intervals as a function of gamma -- #
# --------------------------------------------------------------- #

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
# plot_df_errorbars_true <- data.frame(gamma = gamma_check,
#                                      est_true = rep(estimate_true, length(gamma_check)),
#                                      est_true_lb = rep(estimate_true_lb, length(gamma_check)),
#                                      est_true_ub = rep(estimate_true_ub, length(gamma_check)))

legend_values <- c('avg_xi_hat' = 'palevioletred3',
                   #'est_true' = 'cadetblue',
                   'bounds' = 'black')
legend_values_linetype <- c('avg_xi_hat' = 'solid')
                            #'est_true' = 'solid')
legend_labels <- c('avg_xi_hat' = TeX('$\\hat{\\xi}(A_{\\gamma}^{(tr)})$'))
                   #'est_true' = TeX('$\\hat{\\theta}(A)$ using $Z^{(true)}$'))

dolphins_CI_gamma <- ggplot() +
  geom_line(data = plot_df_midpoint, aes(x = gamma, y = 0), linetype = 'dashed', alpha = 0.6, linewidth = thin_linewidth) +
  geom_line(data = plot_df_midpoint, aes(x = gamma, y = avg_xi_hat, color = 'avg_xi_hat', linetype = 'avg_xi_hat'), linewidth = thick_linewidth) +
  #geom_line(data = plot_df_midpoint, aes(x = gamma, y = estimate_true, color = 'est_true', linetype = 'est_true'), linewidth = thick_linewidth) +
  #geom_linerange(data = plot_df_errorbars_true, aes(x = gamma, y = est_true, ymin = est_true_lb, ymax = est_true_ub, linetype = 'est_true'), color = 'black', linewidth = thin_linewidth, show.legend = FALSE) +
  geom_linerange(data = plot_df_errorbars, aes(x = gamma, y = avg_xi_hat, ymin = lowpoint, ymax = highpoint, linetype = 'est_true'), color = 'black', linewidth = thin_linewidth, show.legend = FALSE) +
  xlab(TeX('$\\gamma$')) + ylab('') +
  labs(color = 'Legend') +
  scale_color_manual(
    values = legend_values,
    labels = legend_labels
  ) +
  # scale_linetype_manual(
  #   values = legend_values_linetype,
  #   labels = legend_labels
  # ) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size)) +
  guides(linetype = 'none')
dolphins_CI_gamma
ggsave('figures/dolphins_CI_gamma.pdf', plot = dolphins_CI_gamma,
       device = 'pdf', width = 14, height = 8, units = 'cm')
