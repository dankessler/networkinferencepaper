library("manynet") # Where dolphins dataset lives
library("networkinference")
library("latex2exp")
library("tidyverse")

ggplot2::theme_set(theme_minimal())

# ====================================
# == Load dolphins adjacency matrix ==
# ====================================

data("ison_dolphins")
dolphins_igraph <- as_igraph(ison_dolphins)

# Get adjacency matrix
dolphins <- igraph::as_adjacency_matrix(dolphins_igraph, sparse = FALSE)
n <- NROW(dolphins)
K <- 2 # Number of estimated communities

# ========================================================== #
# == Conduct inference for selected target, varying gamma == #
# ========================================================== #

gamma_check <- seq(0.008, 0.50, length.out = 25)
num_sim_per_gamma <- 500
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
    if (K == 2) {
      u <- c(1, 0, -2, 1)
    } else if (K == 3) {
      u <- c(1, 0, 0, -1, 1, 0, -1, -1, 1)
    }
    u <- u / sqrt(sum(u^2))

    # Do clustering on the full version of the matrix
    z_hat_full_initial <- nett::spec_clust(dolphins, K = K)
    z_hat_full <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat_full[, i] <- 1 * (z_hat_full_initial == i)
    }
    n_hat_full <- apply(z_hat_full, 2, sum)

    # Split matrix
    dolphins_splitting <- networkinference::split_network(A = dolphins,
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

    # Check agreement of clustering on full and fissioned data
    rand_results_fission_full[gamma_index, rep] <- mclust::adjustedRandIndex(z_hat_full_initial, z_hat_initial)

    # Conduct inference for selected target
    dolphins_inference <- networkinference::infer_network(Ate = dolphins_te,
                                                          u = u,
                                                          communities = z_hat_initial,
                                                          distribution = "bernoulli",
                                                          gamma = gamma,
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

# ==================== #
# == Create figures == #
# ==================== #

axis.text.x.size <- 15
axis.text.y.size <- 15
axis.title.x.size <- 21
axis.title.y.size <- 16
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
plot_df <- data.frame(gamma = gamma_check,
                      avg_rand_fission_full = avg_rand_results_fission_full)

dolphins_randindex_gamma <- ggplot(plot_df) +
  geom_line(aes(x = gamma, y = avg_rand_fission_full), alpha = 0.9, linewidth = thick_linewidth) +
  xlab(TeX('$\\gamma$')) + ylab('Adjusted RAND Index') +
  labs(color = 'Legend') +
  ylim(c(-0.01, 1)) +
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
axis.title.y.size <- 18
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

legend_values <- c('avg_xi_hat' = 'palevioletred3',
                   'bounds' = 'black')
legend_values_linetype <- c('avg_xi_hat' = 'solid')
legend_labels <- c('avg_xi_hat' = TeX('$\\hat{\\xi}(A_{\\gamma}^{(tr)})$'))

dolphins_CI_gamma <- ggplot() +
  geom_line(data = plot_df_midpoint, aes(x = gamma, y = 0), color = 'darkgray', linetype = 'dashed', alpha = 0.8, linewidth = thin_linewidth) +
  geom_line(data = plot_df_midpoint, aes(x = gamma, y = avg_xi_hat), color = 'palevioletred3', linetype = 'solid', linewidth = thick_linewidth) +
  geom_linerange(data = plot_df_errorbars, aes(x = gamma, y = avg_xi_hat, ymin = lowpoint, ymax = highpoint), linetype = 'dashed', color = 'black', linewidth = thin_linewidth, show.legend = FALSE) +
  xlab(TeX('$\\gamma$')) + ylab(TeX('90% CI for $\\xi(A^{(tr)}_\\gamma)$')) +
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
