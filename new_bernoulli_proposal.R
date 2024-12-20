library('networkinference')
library('latex2exp')
library('tidyverse')

ggplot2::theme_set(theme_minimal())

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

# ===========================
# == Simulation parameters ==
# ===========================

n <- 150                                         # Network nodes
K <- 2                                          # Number of communities
num_sim_per_gamma <- 90                        # Repetitions of the simulation
gamma_check <- seq(0.01, 0.50, length.out = 15) # What gamma to check
alpha <- 0.10                                   # Confidence level

# Marginal mean matrix
M <- matrix(rep(0.5, n^2), nrow = n)
M <- matrix(runif(n^2, min = 0.2, max = 0.8), nrow = n)

# The linear combination vector we test
u <- c(1, -1, -1, 1)
u <- u / sqrt(sum(u^2))

# ============================================
# == Check through multiple values of gamma ==
# ============================================

# Where results are stored
varphi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
theta_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
theta_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
theta_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
theta_covered_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

xi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
xi_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
xi_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
xi_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
xi_covered_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

varphi_covered_by_xi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
varphi_covered_by_theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

Capricorn_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

Aquarius1_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
Aquarius2_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
Aquarius3_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
Aquarius4_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
Aquarius_est_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
Aquarius_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))
Aquarius_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma))

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

    # # Clustering on all of A
    # z_hat_full_initial <- nett::spec_clust(A, K = K)
    # z_hat_full <- matrix(rep(NA, n*K), nrow = n)
    # for (i in 1:K) {
    #   z_hat_full[, i] <- 1 * (z_hat_full_initial == i)
    # }
    # n_hat_full <- apply(z_hat_full, 2, sum)

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

    # ============
    # == Varphi ==
    # ============

    NN_inv <- diag(1 / diag(t(z_hat) %*% z_hat))
    C <- NN_inv %*% t(z_hat) %*% M %*% z_hat %*% NN_inv
    varphi <- t(u) %*% as.vector(C)

    # ===========
    # == Theta ==
    # ===========

    # Calculate our estimator and its target
    Cmask <- (gamma / (1-gamma))^(2*A_tr - 1)
    Tmat <- M / (M + (1-M) * Cmask)
    B <- NN_inv %*% t(z_hat) %*% Tmat %*% z_hat %*% NN_inv
    Bhat <- NN_inv %*% t(z_hat) %*% A %*% z_hat %*% NN_inv
    theta_hat <- t(u) %*% as.vector(Bhat)
    theta <- t(u) %*% as.vector(B)

    # Variance
    Delta_theta <- (Bhat * (1 - Bhat)) / (t(z_hat) %*% matrix(1, nrow = n, ncol = n) %*% z_hat)
    Sigma_theta <- t(u) %*% diag(as.vector(Delta_theta)) %*% u

    theta_moe <- qnorm(1 - alpha / 2) * sqrt(Sigma_theta)

    # ========
    # == Xi ==
    # ========

    Phi0_hat <- (t(z_hat) %*% (A * (1 - A_tr)) %*% z_hat) /
      (t(z_hat) %*% (1 - A_tr) %*% z_hat)
    Phi1_hat <- (t(z_hat) %*% (A * A_tr) %*% z_hat) /
      (t(z_hat) %*% A_tr %*% z_hat)
    Phi0_star <- (t(z_hat) %*% (Tmat * (1 - A_tr)) %*% z_hat) /
      (t(z_hat) %*% (1 - A_tr) %*% z_hat)
    Phi1_star <- (t(z_hat) %*% (Tmat * A_tr) %*% z_hat) /
      (t(z_hat) %*% A_tr %*% z_hat)

    # Transform back to original scale
    xi0_hat_components <- h0_inv(Phi0_hat, gamma)
    xi1_hat_components <- h1_inv(Phi1_hat, gamma)
    xi0_components <- h0_inv(Phi0_star, gamma)
    xi1_components <- h1_inv(Phi1_star, gamma)

    xi0 <- t(u) %*% as.vector(xi0_components)
    xi1 <- t(u) %*% as.vector(xi1_components)
    xi0_hat <- t(u) %*% as.vector(xi0_hat_components)
    xi1_hat <- t(u) %*% as.vector(xi1_hat_components)

    # Variance
    Delta0_hat <- (Phi0_hat * (1 - Phi0_hat)) / (t(z_hat) %*% (1 - A_tr) %*% z_hat)
    Delta1_hat <- (Phi1_hat * (1 - Phi1_hat)) / (t(z_hat) %*% A_tr %*% z_hat)
    Sigma0_hat_components <- Delta0_hat * (h0_inv_deriv(Phi0_hat, gamma))^2
    Sigma1_hat_components <- Delta1_hat * (h1_inv_deriv(Phi1_hat, gamma))^2

    Sigma0 <- t(u) %*% diag(as.vector(Sigma0_hat_components)) %*% u
    Sigma1 <- t(u) %*% diag(as.vector(Sigma1_hat_components)) %*% u

    # Optional weighted version
    w0 <- sum(A_tr == 0) / n^2
    w1 <- 1 - w0

    xi <- w0*xi0 + w1*xi1
    xi_hat <- w0*xi0_hat + w1*xi1_hat
    Sigma <- w0^2 * Sigma0 + w1^2 * Sigma1
    xi_moe <- qnorm(1 - alpha / 2) * sqrt(Sigma)

    # ======================================
    # == Aquarius (modification of theta) ==
    # ======================================

    c0 <- log(gamma / (1-gamma))
    taylor_adjustment <- matrix(0, nrow = K, ncol = K)
    drastic_taylor_adjustment <- matrix(0, nrow = K, ncol = K)
    mean_taylor_adjustment <- matrix(0, nrow = K, ncol = K)
    estimate_mean_taylor_adjustment <- matrix(0, nrow = K, ncol = K)
    for (k in 1:K) {
      for (l in 1:K) {
        k_nodes <- which(z_hat_initial == k)
        l_nodes <- which(z_hat_initial == l)

        # Appropriate subseting
        A_tr_subs <- A_tr[k_nodes, l_nodes]
        M_subs <- M[k_nodes, l_nodes]
        mean_M_subs <- mean(M[k_nodes, l_nodes])
        est_M_subs <- mean(A[k_nodes, l_nodes])

        # Cardinality of Ikl, Ikl0, and Ikl1
        Ikl_card <- length(k_nodes) * length(l_nodes)
        Ikl0_card <- sum(A_tr_subs)
        Ikl1_card <- sum(1 - A_tr_subs)

        taylor_adjustment[k, l] <-
          c0 * (sum((1-A_tr_subs) * M_subs * (1 - M_subs)) - sum(A_tr_subs * M_subs * (1 - M_subs))) / Ikl_card
        mean_taylor_adjustment[k, l] <-
          c0 * (sum((1-A_tr_subs) * mean_M_subs * (1 - mean_M_subs)) - sum(A_tr_subs * mean_M_subs * (1 - mean_M_subs))) / Ikl_card
        estimate_mean_taylor_adjustment[k, l] <-
          c0 * (sum((1-A_tr_subs) * est_M_subs * (1 - est_M_subs)) - sum(A_tr_subs * est_M_subs * (1 - est_M_subs))) / Ikl_card
        drastic_taylor_adjustment[k, l] <-
          c0 * (sum((1-A_tr_subs) * 0.25) - sum(A_tr_subs * 0.25)) / Ikl_card
      }
    }
    Aquarius_est <- t(u) %*% as.vector(Bhat - estimate_mean_taylor_adjustment)
    Aquarius1 <- t(u) %*% as.vector(B - estimate_mean_taylor_adjustment)
    Aquarius2 <- t(u) %*% as.vector(B - drastic_taylor_adjustment)
    Aquarius3 <- t(u) %*% as.vector(B - mean_taylor_adjustment)
    Aquarius4 <- t(u) %*% as.vector(B - taylor_adjustment)

    # ===============
    # == Capricorn ==
    # ===============

    # Calculate offsets in logistic regression
    cij_offsets <- (2*A_tr - 1) * log((1-gamma) / gamma)
    CapricornB <- matrix(0, nrow = K, ncol = K)

    for (k in 1:K) {
      for (l in 1:K) {
        k_nodes <- which(z_hat_initial == k)
        l_nodes <- which(z_hat_initial == l)

        A_subset <- A[k_nodes, l_nodes]
        cij_subset <- cij_offsets[k_nodes, l_nodes]

        kl_glm <- glm(as.vector(A_subset) ~ 1 + offset(as.vector(cij_subset)), family = 'binomial')
        CapricornB[k, l] <- expit(kl_glm$coefficients)

        # Sigma0_hat_components <-

        # Variance
        # Delta_hat[k, l] <- (Bhat[k, l] * (1 - Bhat[k, l])) / (length(k_nodes) * length(l_nodes))
        #Delta2_hat[k, l] <- Delta_hat[k, l] * (h)^2
      }
    }
    Capricorn <- t(u) %*% as.vector(CapricornB)

    # ====================
    # == Record coverge ==
    # ====================

    varphi_record[gamma_index, rep] <- varphi

    theta_record[gamma_index, rep] <- theta
    theta_hat_record[gamma_index, rep] <- theta_hat
    theta_ci_lower_record[gamma_index, rep] <- theta_hat - theta_moe
    theta_ci_upper_record[gamma_index, rep] <- theta_hat + theta_moe
    theta_covered_record[gamma_index, rep] <- (theta_hat - theta_moe <= theta) & (theta <= theta_hat + theta_moe)

    xi_record[gamma_index, rep] <- xi
    xi_hat_record[gamma_index, rep] <- xi_hat
    xi_ci_lower_record[gamma_index, rep] <- xi_hat - xi_moe
    xi_ci_upper_record[gamma_index, rep] <- xi_hat + xi_moe
    xi_covered_record[gamma_index, rep] <- (xi_hat - xi_moe <= xi) & (xi <= xi_hat + xi_moe)

    varphi_covered_by_xi_record[gamma_index, rep] <- (xi_hat - xi_moe <= varphi) & (varphi <= xi_hat + xi_moe)
    varphi_covered_by_theta_record[gamma_index, rep] <- (theta_hat - theta_moe <= varphi) & (varphi <= theta_hat + theta_moe)

    Capricorn_record[gamma_index, rep] <- Capricorn

    Aquarius1_record[gamma_index, rep] <- Aquarius1
    Aquarius2_record[gamma_index, rep] <- Aquarius2
    Aquarius3_record[gamma_index, rep] <- Aquarius3
    Aquarius4_record[gamma_index, rep] <- Aquarius4
    Aquarius_est_record[gamma_index, rep] <- Aquarius
    Aquarius_ci_lower_record[gamma_index, rep] <- Aquarius_est - theta_moe
    Aquarius_ci_upper_record[gamma_index, rep] <- Aquarius_est + theta_moe
  }
}

# Check coverage?
# ---------------
# Xi is covered by itself
apply(xi_covered_record, 1, mean)
# Theta is covered by itself
apply(theta_covered_record, 1, mean)
# Varphi is covered by xi
apply(varphi_covered_by_xi_record, 1, mean)
# Varphi is covered by theta
apply(varphi_covered_by_theta_record, 1, mean)

# ====================
# == Create figures ==
# ====================

# -------------------------------------------------
# -- Confidence intervals as a function of gamma --
# -------------------------------------------------

avg_xi <- apply(xi_record, 1, mean)
avg_xi_hat <- apply(xi_hat_record, 1, mean)
avg_xi_lower <- apply(xi_ci_lower_record, 1, mean)
avg_xi_upper <- apply(xi_ci_upper_record, 1, mean)

avg_varphi <- apply(varphi_record, 1, mean)

avg_theta <- apply(theta_record, 1, mean)
avg_theta_hat <- apply(theta_hat_record, 1, mean)
avg_theta_lower <- apply(theta_ci_lower_record, 1, mean)
avg_theta_upper <- apply(theta_ci_upper_record, 1, mean)

avg_Capricorn <- apply(Capricorn_record, 1, mean)
avg_Aquarius_hat <- apply(Aquarius_est_record, 1, mean)
avg_Aquarius1 <- apply(Aquarius1_record, 1, mean)
avg_Aquarius2 <- apply(Aquarius2_record, 1, mean)
avg_Aquarius3 <- apply(Aquarius3_record, 1, mean)
avg_Aquarius4 <- apply(Aquarius4_record, 1, mean)

plot_df <- data.frame(gamma = gamma_check,
                      avg_xi = avg_xi,
                      avg_xi_hat = avg_xi_hat,
                      avg_xi_lower = avg_xi_lower,
                      avg_xi_upper = avg_xi_upper,
                      avg_varphi = avg_varphi,
                      avg_theta = avg_theta,
                      avg_theta_hat = avg_theta_hat,
                      avg_theta_lower = avg_theta_lower,
                      avg_theta_upper = avg_theta_upper,
                      avg_Capricorn = avg_Capricorn,
                      avg_Aquarius_hat = avg_Aquarius_hat,
                      avg_Aquarius1 = avg_Aquarius1,
                      avg_Aquarius2 = avg_Aquarius2,
                      avg_Aquarius3 = avg_Aquarius3,
                      avg_Aquarius4 = avg_Aquarius4)
plot_df2 <- plot_df %>% pivot_longer(c('avg_xi', 'avg_xi_hat',
                                       'avg_xi_lower', 'avg_xi_upper',
                                       'avg_varphi', 'avg_theta', 'avg_theta_hat',
                                       'avg_theta_lower', 'avg_theta_upper', 'avg_Capricorn',
                                       'avg_Aquarius_hat', 'avg_Aquarius1', 'avg_Aquarius2',
                                       'avg_Aquarius3', 'avg_Aquarius4'))

# Cut out some variables
#subsss <=
plot_df2 <- plot_df2[!(plot_df2$name %in%
                        c('avg_xi_lower',
                          'avg_xi_upper',
                          'avg_theta_lower',
                          'avg_theta_upper',
                          'avg_Aquarius_hat',
                          'avg_Aquarius2')), ]


A_CI_gamma <- ggplot(plot_df2) +
  geom_line(aes(x = gamma, y = value, color = name), linewidth = 0.8, alpha = 0.6) +
  xlab(TeX('$\\gamma$')) + ylab('') + # ylab(TeX('$\\theta(A^{(tr)}_\\gamma)$')) +
  labs(color = 'Legend') +
  # scale_color_manual(
  #   labels = c('avg_xi' = TeX('$\\xi$'))
  # ) +
  scale_color_manual(
    values = c('avg_xi' = 'chartreuse3',
               'avg_xi_hat' = 'chartreuse4',
               'avg_varphi' = 'gold3',
               'avg_theta' = 'coral2',
               'avg_theta_hat' = 'coral4',
               'avg_Capricorn' = 'magenta',
               'avg_theta_lower' = 'darkred',
               'avg_theta_upper' = 'darkred',
               'avg_Aquarius_hat' = 'slateblue',
               'avg_Aquarius1' = 'slateblue1',
               'avg_Aquarius2' = 'slateblue2',
               'avg_Aquarius3' = 'slateblue3',
               'avg_Aquarius4' = 'slateblue4',
               'avg_xi_lower' = 'blue4',
               'avg_xi_upper' = 'blue4'),
    # 'bounds2' = 'darkslategray2')
    labels = c('avg_xi' = TeX('$\\xi(A^{(tr)})$'),
               'avg_xi_hat' = TeX('$\\hat{\\xi}(A^{(tr)})$'),
               'avg_varphi' = TeX('$\\varphi(A^{(tr)})$'),
               'avg_theta' = TeX('$\\theta(A^{(tr)})$'),
               'avg_theta_hat' = TeX('$\\hat{\\theta}(A^{(tr)})$'),
               'avg_Capricorn' = TeX('$\\CapricornHat(A^{(tr)})$'),
               'avg_theta_lower' = TeX('$\\theta$ CI lower bound'),
               'avg_theta_upper' = TeX('$\\theta$ CI upper bound'),
               'avg_Aquarius_hat' = TeX('$AquariusHat(A^{(tr)})$'),
               'avg_Aquarius1' = TeX('$Aquarius1(A^{(tr)})$'),
               'avg_Aquarius2' = TeX('$Aquarius2(A^{(tr)})$'),
               'avg_Aquarius3' = TeX('$Aquarius3(A^{(tr)})$'),
               'avg_Aquarius4' = TeX('$Aquarius4(A^{(tr)})$'),
               'avg_xi_lower' = TeX('$\\xi$ CI lower bound'),
               'avg_xi_upper' = TeX('$\\xi$ CI upper bound'))
    #TeX('CI upper bound'))
  ) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
A_CI_gamma

# Find which one is smallest from varphi in MAE
mean(abs(avg_varphi - avg_theta))
mean(abs(avg_varphi - avg_Aquarius1))
mean(abs(avg_varphi - avg_Aquarius2))
mean(abs(avg_varphi - avg_Aquarius3))
mean(abs(avg_varphi - avg_Aquarius4))
mean(abs(avg_varphi - avg_Capricorn))

#' A_CI_gamma <- ggplot(plot_df) +
#'   geom_line(aes(x = gamma, y = avg_xi, color = 'xi'), linewidth = 1.1) +
#'   geom_line(aes(x = gamma, y = avg_xi_hat, color = 'xi_hat'), linewidth = 1.1) +
#'   #geom_line(aes(x = gamma, y = avg_varphi, color = 'varphi'), linewidth = 1.1) +
#'   #geom_line(aes(x = gamma, y = avg_theta, color = 'theta'), linewidth = 1.1) +
#'   # geom_line(aes(x = gamma, y = 0), linetype = 'dashed', alpha = 0.6, linewidth = 0.8) +
#'   geom_line(aes(x = gamma, y = avg_xi_lower, color = 'bounds1'), linewidth = 1.1) +
#'   # geom_line(aes(x = gamma, y = avg_xi_upper, color = 'bounds2'), linewidth = 1.1) +
#'   xlab(TeX('$\\gamma$')) + ylab('') + # ylab(TeX('$\\theta(A^{(tr)}_\\gamma)$')) +
#'   labs(color = 'Legend') +
#'   scale_color_manual(
#'     values = c('xi' = 'firebrick3',
#'                'xi_hat' = 'darkorange',
#'                #'varphi' = 'blue',
#'                #'theta' = 'darkgreen',
#'                'bounds1' = 'darkslategray4'),
#'                # 'bounds2' = 'darkslategray2')
#'     labels = c(TeX('$\\xi(A^{(tr)})$'),
#'                TeX('$\\hat{\\xi}(A^{(tr)})$'),
#'                #TeX('$\\varphi(A^{(tr)})$'),
#'                #TeX('$\\theta(A^{(tr)})$'),
#'                TeX('CI lower bound'))
#'                #TeX('CI upper bound'))
#'   ) +
#'   theme(aspect.ratio = 1,
#'         axis.text.x = element_text(size = 16),
#'         axis.text.y = element_text(size = 16),
#'         axis.title.x = element_text(size = 25),
#'         axis.title.y = element_text(size = 16),
#'         legend.title = element_text(size = 18),
#'         legend.text = element_text(size = 15))
#' A_CI_gamma
# ggsave('figures/bernoulli_drawbacks/fig1.pdf', plot = A_CI_gamma,
#        device = 'pdf', width = 6, height = 4.5)
