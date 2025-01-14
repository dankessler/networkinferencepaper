# =====================================================================
# == Creates figure for demonstrating different targets of inference ==
# =====================================================================

library('tidyverse')
library('latex2exp')
library('patchwork')
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

n <- 150                                          # Network nodes
K <- 2                                            # Number of communities
num_sim_per_gamma_epsilon <- 300                  # Repetitions of the simulation
num_epsilon_gamma_check <- 20
gamma_check <- seq(0.01, 0.50, length.out = num_epsilon_gamma_check)
epsilon_check <- seq(0.01, 0.99, length.out = num_epsilon_gamma_check)
alpha <- 0.10                                     # Confidence level

# The linear combination vector we test
u <- c(1, -1, -1, 1)
u <- u / sqrt(sum(u^2))

# Flat mean matrix of 0.5
M_flat <- matrix(rep(0.5, n^2), nrow = n)
M_flat_poisson <- matrix(rep(2, n^2), nrow = n)

# ================
# == Simulation ==
# ================

# Where results are stored
# ------------------------
# Simulation 1 - Poisson
sim1_poisson_theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_theta_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_varphi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_varphi_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_theta_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_theta_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_varphi_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_varphi_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim1_poisson_theta_coverage <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))

# Simulation 2 - Bernoulli (flat M, "bad target" varphi)
sim2_theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim2_varphi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim2_varphi_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim2_varphi_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim2_varphi_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim2_varphi_coverage <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))

# Simulation 3 - Bernoulli (flat M, better target xi)
sim3_theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim3_xi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim3_xi_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim3_xi_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim3_xi_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim3_xi_coverage <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))

# Simulation 4 - Bernoulli (varying M, better target xi)
sim4_theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim4_xi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim4_xi_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim4_xi_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim4_xi_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim4_xi_coverage <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))

# ------------------------------
# -- Sanity check simulations --
# ------------------------------

# Simulation 5 - Bernoulli (community-based homogeneous M, better target xi)
sim5_theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim5_xi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim5_xi_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim5_xi_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim5_xi_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim5_xi_coverage <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))

# Simulation 5 - Bernoulli (community-based heterogeneous M, better target xi)
sim6_theta_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim6_xi_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim6_xi_hat_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim6_xi_ci_lower_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim6_xi_ci_upper_record <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))
sim6_xi_coverage <- array(0, dim = c(length(gamma_check), num_sim_per_gamma_epsilon))

for (gamma_epsilon_index in 1:num_epsilon_gamma_check) {
  gamma <- gamma_check[gamma_epsilon_index]
  eps <- epsilon_check[gamma_epsilon_index]

  for (rep in 1:num_sim_per_gamma_epsilon) {

    # Simulation 1 drawing
    A_pois <- matrix(rpois(n^2, lambda = as.vector(M_flat_poisson)), nrow = n)
    # Simulations 2, 3 drawing
    A_flat <- matrix(rbinom(n^2, size = 1, prob = as.vector(M_flat)), nrow = n)

    # Simulation 4 drawing
    M_vary_temp <- matrix(rbinom(n = n^2, size = 1, prob = 0.5), nrow = n)
    M_vary <- matrix(0.3, nrow = n, ncol = n) * M_vary_temp + matrix(0.7, nrow = n, ncol = n) * (1 - M_vary_temp)
    A_vary <- matrix(rbinom(n^2, size = 1, prob = as.vector(M_vary)), nrow = n)

    # Simulation 5 and 6 drawing
    # --------------------------
    # With communities appearing, homogeneous/heterogeneous within
    com_hom_means <- matrix(c(0.7, 0.3, 0.3, 0.7), nrow = 2, ncol = 2)
    com_het_means_1 <- matrix(c(0.8, 0.4, 0.4, 0.8), nrow = 2, ncol = 2)
    com_het_means_2 <- matrix(c(0.6, 0.2, 0.2, 0.6), nrow = 2, ncol = 2)
    M_comms_hom <- matrix(rep(NA, n^2), nrow = n)
    M_comms_het <- matrix(rep(NA, n^2), nrow = n)
    for (k in 1:2) {
      for (l in 1:2) {
        # Homogeneous within communities
        M_comms_hom[((k-1)*(n/2) + 1):(k*(n/2)),
                    ((l-1)*(n/2) + 1):(l*(n/2))] <- com_hom_means[k, l]

        # Heterogeneous within communities
        temp <- matrix(rbinom((n/2)^2, size = 1, prob = 0.5), nrow = n/2)
        M_comms_het[((k-1)*(n/2) + 1):(k*(n/2)),
                    ((l-1)*(n/2) + 1):(l*(n/2))] <- temp * com_het_means_1[k, l] + (1 - temp) * com_het_means_2[k, l]
      }
    }
    # Draw A from M
    sim5_A <- matrix(rbinom(n^2, size = 1, prob = as.vector(M_comms_hom)), nrow = n)
    sim6_A <- matrix(rbinom(n^2, size = 1, prob = as.vector(M_comms_het)), nrow = n)

    # ------------------------
    # -- Thinning / Fission --
    # ------------------------

    A_pois_tr <- matrix(rbinom(n = n^2, size = as.vector(A_pois), prob = eps), nrow = n)
    A_pois_te <- A_pois - A_pois_tr

    W <- matrix(rbinom(n^2, size = 1, prob = gamma), nrow = n)
    A_flat_tr <- A_flat * (1 - W) + (1 - A_flat) * W
    A_vary_tr <- A_vary * (1 - W) + (1 - A_vary) * W
    sim5_A_tr <- sim5_A * (1 - W) + (1 - sim5_A) * W
    sim6_A_tr <- sim6_A * (1 - W) + (1 - sim6_A) * W

    # ----------------
    # -- Clustering --
    # ----------------

    # Poisson
    z_hat_pois_initial <- nett::spec_clust(A_pois_tr, K = K)
    z_hat_pois <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat_pois[, i] <- 1 * (z_hat_pois_initial == i)
    }
    n_hat_pois <- apply(z_hat_pois, 2, sum)

    # Bernoulli
    z_hat_flat_initial <- nett::spec_clust(A_flat_tr, K = K)
    z_hat_flat <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat_flat[, i] <- 1 * (z_hat_flat_initial == i)
    }
    n_hat_flat <- apply(z_hat_flat, 2, sum)

    z_hat_vary_initial <- nett::spec_clust(A_vary_tr, K = K)
    z_hat_vary <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      z_hat_vary[, i] <- 1 * (z_hat_vary_initial == i)
    }
    n_hat_vary <- apply(z_hat_vary, 2, sum)

    # Simulation 5 and 6
    sim5_z_hat_initial <- nett::spec_clust(sim5_A_tr, K = K)
    sim5_z_hat <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      sim5_z_hat[, i] <- 1 * (sim5_z_hat_initial == i)
    }
    sim5_n_hat <- apply(sim5_z_hat, 2, sum)

    sim6_z_hat_initial <- nett::spec_clust(sim6_A_tr, K = K)
    sim6_z_hat <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      sim6_z_hat[, i] <- 1 * (sim6_z_hat_initial == i)
    }
    sim6_n_hat <- apply(sim6_z_hat, 2, sum)

    # Some useful matrices
    NN_inv_flat <- diag(1 / diag(t(z_hat_flat) %*% z_hat_flat))
    NN_inv_vary <- diag(1 / diag(t(z_hat_vary) %*% z_hat_vary))
    NN_inv_pois <- diag(1 / diag(t(z_hat_pois) %*% z_hat_pois))
    sim5_NN_inv <- diag(1 / diag(t(sim5_z_hat) %*% sim5_z_hat))
    sim6_NN_inv <- diag(1 / diag(t(sim6_z_hat) %*% sim6_z_hat))

    comm_pair_sample_size_pois <- t(z_hat_pois) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_pois

    comm_pair_sample_size_flat <- t(z_hat_flat) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_flat
    comm_pair_sample_size_flat_ones <- t(z_hat_flat) %*% A_flat_tr %*% z_hat_flat
    comm_pair_sample_size_flat_zeros <- t(z_hat_flat) %*% (1 - A_flat_tr) %*% z_hat_flat

    comm_pair_sample_size_vary <- t(z_hat_vary) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_vary
    comm_pair_sample_size_vary_ones <- t(z_hat_vary) %*% A_vary_tr %*% z_hat_vary
    comm_pair_sample_size_vary_zeros <- t(z_hat_vary) %*% (1 - A_vary_tr) %*% z_hat_vary

    sim5_comm_pair_sample_size <- t(sim5_z_hat) %*% matrix(1, nrow = n, ncol = n) %*% sim5_z_hat
    sim5_comm_pair_sample_size_ones <- t(sim5_z_hat) %*% sim5_A_tr %*% sim5_z_hat
    sim5_comm_pair_sample_size_zeros <- t(sim5_z_hat) %*% (1 - sim5_A_tr) %*% sim5_z_hat

    sim6_comm_pair_sample_size <- t(sim6_z_hat) %*% matrix(1, nrow = n, ncol = n) %*% sim6_z_hat
    sim6_comm_pair_sample_size_ones <- t(sim6_z_hat) %*% sim6_A_tr %*% sim6_z_hat
    sim6_comm_pair_sample_size_zeros <- t(sim6_z_hat) %*% (1 - sim6_A_tr) %*% sim6_z_hat

    # ------------------
    # -- Do inference --
    # ------------------

    # ------------------
    # -- Simulation 1 --
    # ------------------

    sim1_mean_matrix <- NN_inv_pois %*% t(z_hat_pois) %*% M_flat_poisson %*% z_hat_pois %*% NN_inv_pois
    sim1_est_matrix <- NN_inv_pois %*% t(z_hat_pois) %*% A_pois_te %*% z_hat_pois %*% NN_inv_pois
    sim1_var_matrix <- (t(z_hat_pois) %*% ((1-eps) * M_flat_poisson) %*% z_hat_pois) / comm_pair_sample_size_pois^2

    sim1_mean <- t(u) %*% as.vector(sim1_mean_matrix)
    sim1_est <- (1-eps)^(-1) * t(u) %*% as.vector(sim1_est_matrix)
    sim1_moe <- qnorm(1 - alpha / 2) * (1-eps)^(-1) * sqrt(t(u) %*% diag(as.vector(sim1_var_matrix)) %*% u)

    sim1_poisson_theta_record[gamma_epsilon_index, rep] <- sim1_mean
    sim1_poisson_theta_hat_record[gamma_epsilon_index, rep] <- sim1_est
    sim1_poisson_theta_ci_lower_record[gamma_epsilon_index, rep] <- sim1_est - sim1_moe
    sim1_poisson_theta_ci_upper_record[gamma_epsilon_index, rep] <- sim1_est + sim1_moe
    sim1_poisson_theta_coverage[gamma_epsilon_index, rep] <- (sim1_est - sim1_moe <= sim1_mean) & (sim1_mean <= sim1_est + sim1_moe)

    sim1_poisson_varphi_record[gamma_epsilon_index, rep] <- (1-eps) * sim1_mean
    sim1_poisson_varphi_hat_record[gamma_epsilon_index, rep] <- (1-eps) * sim1_est
    sim1_poisson_varphi_ci_lower_record[gamma_epsilon_index, rep] <- (1-eps) * (sim1_est - sim1_moe)
    sim1_poisson_varphi_ci_upper_record[gamma_epsilon_index, rep] <- (1-eps) * (sim1_est + sim1_moe)

    # ------------------
    # -- Simulation 2 --
    # ------------------

    # Calculate T matrix
    Cmask_flat <- (gamma / (1-gamma))^(2*A_flat_tr - 1)
    T_flat <- M_flat / (M_flat + (1-M_flat) * Cmask_flat)

    sim2_marginal_mean_matrix <- NN_inv_flat %*% t(z_hat_flat) %*% M_flat %*% z_hat_flat %*% NN_inv_flat
    sim2_mean_matrix <- NN_inv_flat %*% t(z_hat_flat) %*% T_flat %*% z_hat_flat %*% NN_inv_flat
    sim2_est_matrix <- NN_inv_flat %*% t(z_hat_flat) %*% A_flat %*% z_hat_flat %*% NN_inv_flat
    sim2_var_matrix <- (t(z_hat_flat) %*% (T_flat * (1 - T_flat)) %*% z_hat_flat) / comm_pair_sample_size_flat^2

    sim2_marginal_mean <- t(u) %*% as.vector(sim2_marginal_mean_matrix)
    sim2_mean <- t(u) %*% as.vector(sim2_mean_matrix)
    sim2_est <- t(u) %*% as.vector(sim2_est_matrix)
    sim2_moe <- qnorm(1 - alpha / 2) * sqrt(t(u) %*% diag(as.vector(sim2_var_matrix)) %*% u)

    sim2_theta_record[gamma_epsilon_index, rep] <- sim2_marginal_mean
    sim2_varphi_record[gamma_epsilon_index, rep] <- sim2_mean
    sim2_varphi_hat_record[gamma_epsilon_index, rep] <- sim2_est
    sim2_varphi_ci_lower_record[gamma_epsilon_index, rep] <- sim2_est - sim2_moe
    sim2_varphi_ci_upper_record[gamma_epsilon_index, rep] <- sim2_est + sim2_moe
    sim2_varphi_coverage[gamma_epsilon_index, rep] <- (sim2_est - sim2_moe <= sim2_mean) & (sim2_mean <= sim2_est + sim2_moe)

    # ------------------
    # -- Simulation 3 --
    # ------------------

    # Calculate T matrix
    Cmask_flat <- (gamma / (1-gamma))^(2*A_flat_tr - 1)
    T_flat <- M_flat / (M_flat + (1-M_flat) * Cmask_flat)

    Lambda0_hat_flat <- (t(z_hat_flat) %*% (A_flat * (1 - A_flat_tr)) %*% z_hat_flat) /
      comm_pair_sample_size_flat_zeros
    Lambda1_hat_flat <- (t(z_hat_flat) %*% (A_flat * A_flat_tr) %*% z_hat_flat) /
      comm_pair_sample_size_flat_ones
    Lambda0_star_flat <- (t(z_hat_flat) %*% (T_flat * (1 - A_flat_tr)) %*% z_hat_flat) /
      comm_pair_sample_size_flat_zeros
    Lambda1_star_flat <- (t(z_hat_flat) %*% (T_flat * A_flat_tr) %*% z_hat_flat) /
      comm_pair_sample_size_flat_ones

    Lambda0_var_flat <- (t(z_hat_flat) %*% (T_flat * (1 - T_flat) * (1 - A_flat_tr)) %*% z_hat_flat) /
      comm_pair_sample_size_flat_zeros^2
    Lambda1_var_flat <- (t(z_hat_flat) %*% (T_flat * (1 - T_flat) * A_flat_tr) %*% z_hat_flat) /
      comm_pair_sample_size_flat_ones^2

    # Weighting of 0s and 1s
    ones_weighting_flat <- (t(z_hat_flat) %*% A_flat_tr %*% z_hat_flat) /
      (t(z_hat_flat) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_flat)
    zeros_weighting_flat <- (t(z_hat_flat) %*% (1 - A_flat_tr) %*% z_hat_flat) /
      (t(z_hat_flat) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_flat)

    Phi_matrix <- zeros_weighting_flat * h0_inv(Lambda0_star_flat, gamma) + ones_weighting_flat * h1_inv(Lambda1_star_flat, gamma)
    Phi_hat_matrix <- zeros_weighting_flat * h0_inv(Lambda0_hat_flat, gamma) + ones_weighting_flat * h1_inv(Lambda1_hat_flat, gamma)
    Phi_var_matrix <- zeros_weighting_flat^2 * Lambda0_var_flat * (h0_inv_deriv(Lambda0_star_flat, gamma))^2 +
      ones_weighting_flat^2 * Lambda1_var_flat * (h1_inv_deriv(Lambda1_star_flat, gamma))^2

    xi <- t(u) %*% as.vector(Phi_matrix)
    xi_hat <- t(u) %*% as.vector(Phi_hat_matrix)
    xi_moe <- qnorm(1 - alpha/2) * sqrt(t(u) %*% diag(as.vector(Phi_var_matrix)) %*% u)

    sim3_theta_record[gamma_epsilon_index, rep] <- t(u) %*% as.vector(NN_inv_flat %*% t(z_hat_flat) %*% M_flat %*% z_hat_flat %*% NN_inv_flat)
    sim3_xi_record[gamma_epsilon_index, rep] <- xi
    sim3_xi_hat_record[gamma_epsilon_index, rep] <- xi_hat
    sim3_xi_ci_lower_record[gamma_epsilon_index, rep] <- xi_hat - xi_moe
    sim3_xi_ci_upper_record[gamma_epsilon_index, rep] <- xi_hat + xi_moe
    sim3_xi_coverage[gamma_epsilon_index, rep] <- (xi_hat - xi_moe <= xi) & (xi <= xi_hat + xi_moe)

    # ------------------
    # -- Simulation 4 --
    # ------------------

    # Calculate T matrix
    Cmask_vary <- (gamma / (1-gamma))^(2*A_vary_tr - 1)
    T_vary <- M_vary / (M_vary + (1-M_vary) * Cmask_vary)

    Lambda0_hat_vary <- (t(z_hat_vary) %*% (A_vary * (1 - A_vary_tr)) %*% z_hat_vary) /
      comm_pair_sample_size_vary_zeros
    Lambda1_hat_vary <- (t(z_hat_vary) %*% (A_vary * A_vary_tr) %*% z_hat_vary) /
      comm_pair_sample_size_vary_ones
    Lambda0_star_vary <- (t(z_hat_vary) %*% (T_vary * (1 - A_vary_tr)) %*% z_hat_vary) /
      comm_pair_sample_size_vary_zeros
    Lambda1_star_vary <- (t(z_hat_vary) %*% (T_vary * A_vary_tr) %*% z_hat_vary) /
      comm_pair_sample_size_vary_ones

    Lambda0_var_vary <- (t(z_hat_vary) %*% (T_vary * (1 - T_vary) * (1 - A_vary_tr)) %*% z_hat_vary) /
      comm_pair_sample_size_vary_zeros^2
    Lambda1_var_vary <- (t(z_hat_vary) %*% (T_vary * (1 - T_vary) * A_vary_tr) %*% z_hat_vary) /
      comm_pair_sample_size_vary_ones^2

    # Weighting of 0s and 1s
    ones_weighting_vary <- (t(z_hat_vary) %*% A_vary_tr %*% z_hat_vary) /
      (t(z_hat_vary) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_vary)
    zeros_weighting_vary <- (t(z_hat_vary) %*% (1 - A_vary_tr) %*% z_hat_vary) /
      (t(z_hat_vary) %*% matrix(1, nrow = n, ncol = n) %*% z_hat_vary)

    Phi_matrix <- zeros_weighting_vary * h0_inv(Lambda0_star_vary, gamma) + ones_weighting_vary * h1_inv(Lambda1_star_vary, gamma)
    Phi_hat_matrix <- zeros_weighting_vary * h0_inv(Lambda0_hat_vary, gamma) + ones_weighting_vary * h1_inv(Lambda1_hat_vary, gamma)
    Phi_var_matrix <- zeros_weighting_vary^2 * Lambda0_var_vary * (h0_inv_deriv(Lambda0_star_vary, gamma))^2 +
      ones_weighting_vary^2 * Lambda1_var_vary * (h1_inv_deriv(Lambda1_star_vary, gamma))^2

    xi <- t(u) %*% as.vector(Phi_matrix)
    xi_hat <- t(u) %*% as.vector(Phi_hat_matrix)
    xi_moe <- qnorm(1 - alpha/2) * sqrt(t(u) %*% diag(as.vector(Phi_var_matrix)) %*% u)

    sim4_theta_record[gamma_epsilon_index, rep] <- t(u) %*% as.vector(NN_inv_vary %*% t(z_hat_vary) %*% M_vary %*% z_hat_vary %*% NN_inv_vary)
    sim4_xi_record[gamma_epsilon_index, rep] <- xi
    sim4_xi_hat_record[gamma_epsilon_index, rep] <- xi_hat
    sim4_xi_ci_lower_record[gamma_epsilon_index, rep] <- xi_hat - xi_moe
    sim4_xi_ci_upper_record[gamma_epsilon_index, rep] <- xi_hat + xi_moe
    sim4_xi_coverage[gamma_epsilon_index, rep] <- (xi_hat - xi_moe <= xi) & (xi <= xi_hat + xi_moe)

    # ------------------
    # -- Simulation 5 --
    # ------------------

    # Calculate T matrix
    sim5_Cmask <- (gamma / (1-gamma))^(2*sim5_A_tr - 1)
    sim5_T <- M_comms_hom / (M_comms_hom + (1-M_comms_hom) * sim5_Cmask)

    sim5_Lambda0_hat <- (t(sim5_z_hat) %*% (sim5_A * (1 - sim5_A_tr)) %*% sim5_z_hat) /
      sim5_comm_pair_sample_size_zeros
    sim5_Lambda1_hat <- (t(sim5_z_hat) %*% (sim5_A * sim5_A_tr) %*% sim5_z_hat) /
      sim5_comm_pair_sample_size_ones
    sim5_Lambda0_star <- (t(sim5_z_hat) %*% (sim5_T * (1 - sim5_A_tr)) %*% sim5_z_hat) /
      sim5_comm_pair_sample_size_zeros
    sim5_Lambda1_star <- (t(sim5_z_hat) %*% (sim5_T * sim5_A_tr) %*% sim5_z_hat) /
      sim5_comm_pair_sample_size_ones

    sim5_Lambda0_var <- (t(sim5_z_hat) %*% (sim5_T * (1 - sim5_T) * (1 - sim5_A_tr)) %*% sim5_z_hat) /
      sim5_comm_pair_sample_size_zeros^2
    sim5_Lambda1_var <- (t(sim5_z_hat) %*% (sim5_T * (1 - sim5_T) * sim5_A_tr) %*% sim5_z_hat) /
      sim5_comm_pair_sample_size_ones^2

    # Weighting of 0s and 1s
    sim5_ones_weighting <- (t(sim5_z_hat) %*% sim5_A_tr %*% sim5_z_hat) /
      (t(sim5_z_hat) %*% matrix(1, nrow = n, ncol = n) %*% sim5_z_hat)
    sim5_zeros_weighting <- (t(sim5_z_hat) %*% (1 - sim5_A_tr) %*% sim5_z_hat) /
      (t(sim5_z_hat) %*% matrix(1, nrow = n, ncol = n) %*% sim5_z_hat)

    sim5_Phi_matrix <- sim5_zeros_weighting * h0_inv(sim5_Lambda0_star, gamma) + sim5_ones_weighting * h1_inv(sim5_Lambda1_star, gamma)
    sim5_Phi_hat_matrix <- sim5_zeros_weighting * h0_inv(sim5_Lambda0_hat, gamma) + sim5_ones_weighting * h1_inv(sim5_Lambda1_hat, gamma)
    sim5_Phi_var_matrix <- sim5_zeros_weighting^2 * sim5_Lambda0_var * (h0_inv_deriv(sim5_Lambda0_star, gamma))^2 +
      sim5_ones_weighting^2 * sim5_Lambda1_var * (h1_inv_deriv(sim5_Lambda1_star, gamma))^2

    sim5_xi <- t(u) %*% as.vector(sim5_Phi_matrix)
    sim5_xi_hat <- t(u) %*% as.vector(sim5_Phi_hat_matrix)
    sim5_xi_moe <- qnorm(1 - alpha/2) * sqrt(t(u) %*% diag(as.vector(sim5_Phi_var_matrix)) %*% u)

    sim5_theta_record[gamma_epsilon_index, rep] <- t(u) %*% as.vector(sim5_NN_inv %*% t(sim5_z_hat) %*% M_comms_hom %*% sim5_z_hat %*% sim5_NN_inv)
    sim5_xi_record[gamma_epsilon_index, rep] <- sim5_xi
    sim5_xi_hat_record[gamma_epsilon_index, rep] <- sim5_xi_hat
    sim5_xi_ci_lower_record[gamma_epsilon_index, rep] <- sim5_xi_hat - sim5_xi_moe
    sim5_xi_ci_upper_record[gamma_epsilon_index, rep] <- sim5_xi_hat + sim5_xi_moe
    sim5_xi_coverage[gamma_epsilon_index, rep] <- (sim5_xi_hat - sim5_xi_moe <= sim5_xi) & (sim5_xi <= sim5_xi_hat + sim5_xi_moe)

    # ------------------
    # -- Simulation 6 --
    # ------------------

    # Calculate T matrix
    sim6_Cmask <- (gamma / (1-gamma))^(2*sim6_A_tr - 1)
    sim6_T <- M_comms_hom / (M_comms_hom + (1-M_comms_hom) * sim6_Cmask)

    sim6_Lambda0_hat <- (t(sim6_z_hat) %*% (sim6_A * (1 - sim6_A_tr)) %*% sim6_z_hat) /
      sim6_comm_pair_sample_size_zeros
    sim6_Lambda1_hat <- (t(sim6_z_hat) %*% (sim6_A * sim6_A_tr) %*% sim6_z_hat) /
      sim6_comm_pair_sample_size_ones
    sim6_Lambda0_star <- (t(sim6_z_hat) %*% (sim6_T * (1 - sim6_A_tr)) %*% sim6_z_hat) /
      sim6_comm_pair_sample_size_zeros
    sim6_Lambda1_star <- (t(sim6_z_hat) %*% (sim6_T * sim6_A_tr) %*% sim6_z_hat) /
      sim6_comm_pair_sample_size_ones

    sim6_Lambda0_var <- (t(sim6_z_hat) %*% (sim6_T * (1 - sim6_T) * (1 - sim6_A_tr)) %*% sim6_z_hat) /
      sim6_comm_pair_sample_size_zeros^2
    sim6_Lambda1_var <- (t(sim6_z_hat) %*% (sim6_T * (1 - sim6_T) * sim6_A_tr) %*% sim6_z_hat) /
      sim6_comm_pair_sample_size_ones^2

    # Weighting of 0s and 1s
    sim6_ones_weighting <- (t(sim6_z_hat) %*% sim6_A_tr %*% sim6_z_hat) /
      (t(sim6_z_hat) %*% matrix(1, nrow = n, ncol = n) %*% sim6_z_hat)
    sim6_zeros_weighting <- (t(sim6_z_hat) %*% (1 - sim6_A_tr) %*% sim6_z_hat) /
      (t(sim6_z_hat) %*% matrix(1, nrow = n, ncol = n) %*% sim6_z_hat)

    sim6_Phi_matrix <- sim6_zeros_weighting * h0_inv(sim6_Lambda0_star, gamma) + sim6_ones_weighting * h1_inv(sim6_Lambda1_star, gamma)
    sim6_Phi_hat_matrix <- sim6_zeros_weighting * h0_inv(sim6_Lambda0_hat, gamma) + sim6_ones_weighting * h1_inv(sim6_Lambda1_hat, gamma)
    sim6_Phi_var_matrix <- sim6_zeros_weighting^2 * sim6_Lambda0_var * (h0_inv_deriv(sim6_Lambda0_star, gamma))^2 +
      sim6_ones_weighting^2 * sim6_Lambda1_var * (h1_inv_deriv(sim6_Lambda1_star, gamma))^2

    sim6_xi <- t(u) %*% as.vector(sim6_Phi_matrix)
    sim6_xi_hat <- t(u) %*% as.vector(sim6_Phi_hat_matrix)
    sim6_xi_moe <- qnorm(1 - alpha/2) * sqrt(t(u) %*% diag(as.vector(sim6_Phi_var_matrix)) %*% u)

    sim6_theta_record[gamma_epsilon_index, rep] <- t(u) %*% as.vector(sim6_NN_inv %*% t(sim6_z_hat) %*% M_comms_hom %*% sim6_z_hat %*% sim6_NN_inv)
    sim6_xi_record[gamma_epsilon_index, rep] <- sim6_xi
    sim6_xi_hat_record[gamma_epsilon_index, rep] <- sim6_xi_hat
    sim6_xi_ci_lower_record[gamma_epsilon_index, rep] <- sim6_xi_hat - sim6_xi_moe
    sim6_xi_ci_upper_record[gamma_epsilon_index, rep] <- sim6_xi_hat + sim6_xi_moe
    sim6_xi_coverage[gamma_epsilon_index, rep] <- (sim6_xi_hat - sim6_xi_moe <= sim6_xi) & (sim6_xi <= sim6_xi_hat + sim6_xi_moe)
  }
}

# Sanity check for coverage
apply(sim1_poisson_theta_coverage, 1, mean)
apply(sim2_varphi_coverage, 1, mean)
apply(sim3_xi_coverage, 1, mean)
apply(sim4_xi_coverage, 1, mean)
apply(sim5_xi_coverage, 1, mean)
apply(sim6_xi_coverage, 1, mean)

# ===========
# == Plots ==
# ===========

# --------------------------------
# -- Standardized plot controls --
# --------------------------------

axis.text.x.size <- 13
axis.text.y.size <- 13
axis.title.x.size <- 21
axis.title.y.size <- 21
legend.title.size <- 18
legend.text.size <- 14.5

legend.linewidth.thick <- 0.9
legend.linewidth.thin <- 0.7

thick_linewidth = 1.0
thin_linewidth = 0.9

line_alpha <- 0.9

# Creating a legend that is shared with all six plots
legend_values <- c('avg_theta' = 'slateblue4',
                   'avg_theta_hat' = 'slateblue2',
                   'avg_varphi' = 'slateblue1',
                   'avg_varphi_hat' = 'slateblue4',
                   'avg_xi' = 'palevioletred4',
                   'avg_xi_hat' = 'palevioletred3',
                   'ci_bounds_theta' = 'slateblue2',
                   'ci_bounds_varphi' = 'slateblue4',
                   'ci_bounds_xi' = 'palevioletred3')
legend_values_linetype <- c('avg_theta' = 'solid',
                            'avg_theta_hat' = 'solid',
                            'avg_varphi' = 'solid',
                            'avg_varphi_hat' = 'solid',
                            'avg_xi' = 'solid',
                            'avg_xi_hat' = 'solid',
                            'ci_bounds_theta' = 'dashed',
                            'ci_bounds_varphi' = 'dashed',
                            'ci_bounds_xi' = 'dashed')
legend_labels <- c('avg_theta' = TeX('$\\theta(A^{(tr)})$'),
                   'avg_theta_hat' = TeX('$\\hat{\\theta}(A^{(tr)})$'),
                   'avg_varphi' = TeX('$\\varphi(A^{(tr)})$'),
                   'avg_varphi_hat' = TeX('$\\hat{\\varphi}(A^{(tr)})$'),
                   'avg_xi' = TeX('$\\xi(A^{(tr)})$'),
                   'avg_xi_hat' = TeX('$\\hat{\\xi}(A^{(tr)})$'),
                   'ci_bounds_theta' = TeX('90% CI for $\\theta(A^{(tr)})$'),
                   'ci_bounds_varphi' = TeX('90% CI for $\\varphi(A^{(tr)})$'),
                   'ci_bounds_xi' = TeX('90% CI for $\\xi(A^{(tr)})$'))

unified_plot_levels <- c('avg_theta', 'avg_theta_hat', 'avg_varphi', 'avg_varphi_hat', 'avg_xi', 'avg_xi_hat', 'ci_bounds_theta', 'ci_bounds_varphi', 'ci_bounds_xi')
shared_color_scale <- scale_color_manual(values = legend_values,
                                         labels = legend_labels,
                                         breaks = unified_plot_levels)
shared_linetype_scale <- scale_linetype_manual(values = legend_values_linetype,
                                               labels = legend_labels,
                                               breaks = unified_plot_levels)

# Creating a legend that separates into the three different rows and only
# has the relevant variables
legend_values_row1 <- c('avg_theta' = 'slateblue4',
                        'avg_theta_hat' = 'slateblue2',
                        'avg_varphi' = 'darkgreen',
                        'avg_varphi_hat' = 'darkseagreen3',
                        'ci_bounds_theta' = 'slateblue2',
                        'ci_bounds_varphi' = 'darkseagreen3')
legend_values_linetype_row1 <- c('avg_theta' = 'solid',
                                 'avg_theta_hat' = 'solid',
                                 'avg_varphi' = 'solid',
                                 'avg_varphi_hat' = 'solid',
                                 'ci_bounds_theta' = 'dashed',
                                 'ci_bounds_varphi' = 'dashed')
legend_labels_row1 <- c('avg_theta' = TeX('$\\theta(A^{(tr)})$'),
                        'avg_theta_hat' = TeX('$\\hat{\\theta}(A^{(tr)})$ (left only)'),
                        'avg_varphi' = TeX('$\\varphi(A^{(tr)})$'),
                        'avg_varphi_hat' = TeX('$\\hat{\\varphi}(A^{(tr)})$'),
                        'ci_bounds_theta' = TeX('90% CI for $\\theta(A^{(tr)})$ (left only)'),
                        'ci_bounds_varphi' = TeX('90% CI for $\\varphi(A^{(tr)})$'))
unified_plot_levels_row1 <- c('avg_theta', 'avg_theta_hat', 'avg_varphi', 'avg_varphi_hat', 'ci_bounds_theta', 'ci_bounds_varphi')
shared_color_scale_row1 <- scale_color_manual(values = legend_values_row1,
                                              labels = legend_labels_row1,
                                              breaks = unified_plot_levels_row1)
shared_linetype_scale_row1 <- scale_linetype_manual(values = legend_values_linetype_row1,
                                                    labels = legend_labels_row1,
                                                    breaks = unified_plot_levels_row1)

legend_values_row2 <- c('avg_theta' = 'slateblue4',
                        'avg_xi' = 'palevioletred4',
                        'avg_xi_hat' = 'palevioletred3',
                        'ci_bounds_xi' = 'palevioletred3')
legend_labels_row2 <- c('avg_theta' = TeX('$\\theta(A^{(tr)})$'),
                        'avg_xi' = TeX('$\\xi(A^{(tr)})$'),
                        'avg_xi_hat' = TeX('$\\hat{\\xi}(A^{(tr)})$'),
                        'ci_bounds_xi' = TeX('90% CI for $\\xi(A^{(tr)})$'))
legend_values_linetype_row2 <- c('avg_theta' = 'solid',
                                 'avg_xi' = 'solid',
                                 'avg_xi_hat' = 'solid',
                                 'ci_bounds_xi' = 'dashed')
unified_plot_levels_row2 <- c('avg_theta', 'avg_xi', 'avg_xi_hat', 'ci_bounds_xi')
shared_color_scale_row2 <- scale_color_manual(values = legend_values_row2,
                                              labels = legend_labels_row2,
                                              breaks = unified_plot_levels_row2)
shared_linetype_scale_row2 <- scale_linetype_manual(values = legend_values_linetype_row2,
                                                    labels = legend_labels_row2,
                                                    breaks = unified_plot_levels_row2)

legend_values_row3 <- c('avg_theta' = 'slateblue4',
                        'avg_xi' = 'palevioletred4',
                        'avg_xi_hat' = 'palevioletred3',
                        'ci_bounds_xi' = 'palevioletred3')
legend_labels_row3 <- c('avg_theta' = TeX('$\\theta(A^{(tr)})$'),
                        'avg_xi' = TeX('$\\xi(A^{(tr)})$'),
                        'avg_xi_hat' = TeX('$\\hat{\\xi}(A^{(tr)})$'),
                        'ci_bounds_xi' = TeX('90% CI for $\\xi(A^{(tr)})$'))
legend_values_linetype_row3 <- c('avg_theta' = 'solid',
                                 'avg_xi' = 'solid',
                                 'avg_xi_hat' = 'solid',
                                 'ci_bounds_xi' = 'dashed')
unified_plot_levels_row3 <- c('avg_theta', 'avg_xi', 'avg_xi_hat', 'ci_bounds_xi')
shared_color_scale_row3 <- scale_color_manual(values = legend_values_row3,
                                              labels = legend_labels_row3,
                                              breaks = unified_plot_levels_row3)
shared_linetype_scale_row3 <- scale_linetype_manual(values = legend_values_linetype_row3,
                                                    labels = legend_labels_row3,
                                                    breaks = unified_plot_levels_row3)

# Averages of all relevant variables
sim1_avg_theta <- apply(sim1_poisson_theta_record, 1, mean)
sim1_avg_theta_hat <- apply(sim1_poisson_theta_hat_record, 1, mean)
sim1_avg_theta_lower <- apply(sim1_poisson_theta_ci_lower_record, 1, mean)
sim1_avg_theta_upper <- apply(sim1_poisson_theta_ci_upper_record, 1, mean)
sim1_avg_varphi <- apply(sim1_poisson_varphi_record, 1, mean)
sim1_avg_varphi_hat <- apply(sim1_poisson_varphi_hat_record, 1, mean)
sim1_avg_varphi_lower <- apply(sim1_poisson_varphi_lower_record, 1, mean)
sim1_avg_varphi_upper <- apply(sim1_poisson_varphi_upper_record, 1, mean)

sim2_avg_varphi <- apply(sim2_varphi_record, 1, mean)
sim2_avg_varphi_hat <- apply(sim2_varphi_hat_record, 1, mean)
sim2_avg_varphi_lower <- apply(sim2_varphi_ci_lower_record, 1, mean)
sim2_avg_varphi_upper <- apply(sim2_varphi_ci_upper_record, 1, mean)
sim2_avg_theta <- apply(sim2_theta_record, 1, mean)

sim3_avg_xi <- apply(sim3_xi_record, 1, mean)
sim3_avg_xi_hat <- apply(sim3_xi_hat_record, 1, mean)
sim3_avg_xi_lower <- apply(sim3_xi_ci_lower_record, 1, mean)
sim3_avg_xi_upper <- apply(sim3_xi_ci_upper_record, 1, mean)
sim3_avg_theta <- apply(sim3_theta_record, 1, mean)

sim4_avg_xi <- apply(sim4_xi_record, 1, mean)
sim4_avg_xi_hat <- apply(sim4_xi_hat_record, 1, mean)
sim4_avg_xi_lower <- apply(sim4_xi_ci_lower_record, 1, mean)
sim4_avg_xi_upper <- apply(sim4_xi_ci_upper_record, 1, mean)
sim4_avg_theta <- apply(sim4_theta_record, 1, mean)

sim5_avg_xi <- apply(sim5_xi_record, 1, mean)
sim5_avg_xi_hat <- apply(sim5_xi_hat_record, 1, mean)
sim5_avg_xi_lower <- apply(sim5_xi_ci_lower_record, 1, mean)
sim5_avg_xi_upper <- apply(sim5_xi_ci_upper_record, 1, mean)
sim5_avg_theta <- apply(sim5_theta_record, 1, mean)

sim6_avg_xi <- apply(sim6_xi_record, 1, mean)
sim6_avg_xi_hat <- apply(sim6_xi_hat_record, 1, mean)
sim6_avg_xi_lower <- apply(sim6_xi_ci_lower_record, 1, mean)
sim6_avg_xi_upper <- apply(sim6_xi_ci_upper_record, 1, mean)
sim6_avg_theta <- apply(sim6_theta_record, 1, mean)

# ------------
# -- Plot 1 --
# ------------
# Poisson var phi
plot_df <- data.frame(eps_gamma = 1 - epsilon_check,
                      avg_theta = sim1_avg_theta,
                      avg_theta_hat = sim1_avg_theta_hat,
                      avg_varphi = sim1_avg_varphi,
                      avg_varphi_hat = sim1_avg_varphi_hat) %>%
  pivot_longer(c('avg_theta', 'avg_theta_hat', 'avg_varphi', 'avg_varphi_hat')) %>%
  complete(name = unified_plot_levels_row1)

plot_df_ci_lower <- data.frame(eps_gamma = 1 - epsilon_check,
                               ci_bounds_theta = sim1_avg_theta_lower) %>%
  pivot_longer(c('ci_bounds_theta')) %>%
  complete(name = unified_plot_levels_row1)
plot_df_ci_upper <- data.frame(eps_gamma = 1 - epsilon_check,
                               ci_bounds_theta = sim1_avg_theta_upper) %>%
  pivot_longer(c('ci_bounds_theta')) %>%
  complete(name = unified_plot_levels_row1)
plot_df_ci_varphi_lower <- data.frame(eps_gamma = 1 - epsilon_check,
                               ci_bounds_varphi = sim1_avg_varphi_lower) %>%
  pivot_longer(c('ci_bounds_theta')) %>%
  complete(name = unified_plot_levels_row1)
plot_df_ci_varphi_upper <- data.frame(eps_gamma = 1 - epsilon_check,
                               ci_bounds_varphi = sim1_avg_varphi_upper) %>%
  pivot_longer(c('ci_bounds_theta')) %>%
  complete(name = unified_plot_levels_row1)

# Give consistent factor levels to each of the plot_dfs
plot_df$name <- factor(plot_df$name, levels = unified_plot_levels_row1)
plot_df_ci_lower$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row1)
plot_df_ci_upper$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row1)

plot1 <- ggplot() +
  geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thick_linewidth, alpha = line_alpha) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name)) +
  # geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  geom_line(data = plot_df_ci_lower, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  geom_line(data = plot_df_ci_upper, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  xlab(TeX('$1 - \\epsilon$')) + ylab('') +
  labs(color = 'Legend') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size))

# Add labels
plot1_row_form <- plot1 + shared_color_scale_row1 + shared_linetype_scale_row1 +
  guides(color = guide_legend(override.aes = list(
    linetype = c('solid', 'solid', 'solid', 'solid', 'dashed', 'dashed'),
    linewidth = c(legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thin, legend.linewidth.thin)
  )),
  linetype = 'none')
# plot1 <- plot1 + shared_color_scale + shared_linetype_scale +
#   guides(color = guide_legend(override.aes = list(
#     linetype = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed'),
#     linewidth = 0.6
#   )),
#   linetype = 'none')

plot1_row_form

# ------------
# -- Plot 2 --
# ------------
# Bernoulli varphi
plot_df <- data.frame(eps_gamma = gamma_check,
                      avg_varphi = sim2_avg_varphi,
                      avg_varphi_hat = sim2_avg_varphi_hat,
                      avg_theta = sim2_avg_theta) %>%
  pivot_longer(c('avg_varphi', 'avg_varphi_hat', 'avg_theta')) %>%
  complete(name = unified_plot_levels_row1)
plot_df_ci_lower <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_varphi = sim2_avg_varphi_lower) %>%
  pivot_longer(c('ci_bounds_varphi')) %>%
  complete(name = unified_plot_levels_row1)
plot_df_ci_upper <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_varphi = sim2_avg_varphi_upper) %>%
  pivot_longer(c('ci_bounds_varphi')) %>%
  complete(name = unified_plot_levels_row1)

# Give consistent factor levels to each of the plot_dfs
plot_df$name <- factor(plot_df$name, levels = unified_plot_levels_row1)
plot_df_ci_lower$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row1)
plot_df_ci_upper$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row1)

plot2 <- ggplot() +
  geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thick_linewidth, alpha = line_alpha) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name)) +
  # geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  geom_line(data = plot_df_ci_lower, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  geom_line(data = plot_df_ci_upper, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  xlab(TeX('$\\gamma$')) + ylab('') +
  labs(color = 'Legend') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size))

# Add labels
plot2_row_form <- plot2 + shared_color_scale_row1 + shared_linetype_scale_row1 +
  guides(color = guide_legend(override.aes = list(
    linetype = c('solid', 'solid', 'solid', 'solid', 'dashed', 'dashed'),
    linewidth = c(legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thin, legend.linewidth.thin)
  )),
  linetype = 'none')
# plot2 <- plot2 + shared_color_scale + shared_linetype_scale +
#   guides(color = guide_legend(override.aes = list(
#     linetype = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed'),
#     linewidth = 0.6
#   )),
#   linetype = 'none')

plot2_row_form

# ------------
# -- Plot 3 --
# ------------
# Bernoulli xi
plot_df <- data.frame(eps_gamma = gamma_check,
                      avg_xi = sim3_avg_xi,
                      avg_xi_hat = sim3_avg_xi_hat,
                      avg_theta = sim3_avg_theta) %>%
  pivot_longer(c('avg_xi', 'avg_xi_hat', 'avg_theta')) %>%
  complete(name = unified_plot_levels_row2)
plot_df_ci_lower <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim3_avg_xi_lower) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row2)
plot_df_ci_upper <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim3_avg_xi_upper) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row2)

# Give consistent factor levels to each of the plot_dfs
plot_df$name <- factor(plot_df$name, levels = unified_plot_levels_row2)
plot_df_ci_lower$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row2)
plot_df_ci_upper$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row2)

plot3 <- ggplot() +
  geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thick_linewidth, alpha = line_alpha) +
  geom_line(data = plot_df_ci_lower, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  geom_line(data = plot_df_ci_upper, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  xlab(TeX('$\\gamma$')) + ylab('') +
  labs(color = 'Legend') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size))

# Add labels
plot3_row_form <- plot3 + shared_color_scale_row2 + shared_linetype_scale_row2 +
  guides(color = guide_legend(override.aes = list(
    linetype = c('solid', 'solid', 'solid', 'dashed'),
    linewidth = c(legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thin)
  )),
  linetype = 'none')
# plot3 <- plot3 + shared_color_scale + shared_linetype_scale +
#   guides(color = guide_legend(override.aes = list(
#     linetype = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed'),
#     linewidth = 0.6
#   )),
#   linetype = 'none')

plot3_row_form

# ------------
# -- Plot 4 --
# ------------
# Bernoulli xi
plot_df <- data.frame(eps_gamma = gamma_check,
                      avg_xi = sim4_avg_xi,
                      avg_xi_hat = sim4_avg_xi_hat,
                      avg_theta = sim4_avg_theta) %>%
  pivot_longer(c('avg_xi', 'avg_xi_hat', 'avg_theta')) %>%
  complete(name = unified_plot_levels_row2)
plot_df_ci_lower <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim4_avg_xi_lower) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row2)
plot_df_ci_upper <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim4_avg_xi_upper) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row2)

# Give consistent factor levels to each of the plot_dfs
plot_df$name <- factor(plot_df$name, levels = unified_plot_levels_row2)
plot_df_ci_lower$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row2)
plot_df_ci_upper$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row2)

plot4 <- ggplot() +
  geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thick_linewidth, alpha = line_alpha) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name)) +
  # geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  geom_line(data = plot_df_ci_lower, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  geom_line(data = plot_df_ci_upper, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  xlab(TeX('$\\gamma$')) + ylab('') +
  labs(color = 'Legend') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size))

# Add labels
plot4_row_form <- plot4 + shared_color_scale_row2 + shared_linetype_scale_row2 +
  guides(color = guide_legend(override.aes = list(
    linetype = c('solid', 'solid', 'solid', 'dashed'),
    linewidth = c(legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thin)
  )),
  linetype = 'none')
# plot4 <- plot4 + shared_color_scale + shared_linetype_scale +
#   guides(color = guide_legend(override.aes = list(
#     linetype = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed'),
#     linewidth = 0.6
#   )),
#   linetype = 'none')

plot4_row_form

# ------------
# -- Plot 5 --
# ------------
# Bernoulli xi
plot_df <- data.frame(eps_gamma = gamma_check,
                      avg_xi = sim5_avg_xi,
                      avg_xi_hat = sim5_avg_xi_hat,
                      avg_theta = sim5_avg_theta) %>%
  pivot_longer(c('avg_xi', 'avg_xi_hat', 'avg_theta')) %>%
  complete(name = unified_plot_levels_row3)
plot_df_ci_lower <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim5_avg_xi_lower) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row3)
plot_df_ci_upper <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim5_avg_xi_upper) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row3)

# Give consistent factor levels to each of the plot_dfs
plot_df$name <- factor(plot_df$name, levels = unified_plot_levels_row3)
plot_df_ci_lower$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row3)
plot_df_ci_upper$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row3)

plot5 <- ggplot() +
  geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thick_linewidth, alpha = line_alpha) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name)) +
  # geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  geom_line(data = plot_df_ci_lower, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  geom_line(data = plot_df_ci_upper, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  xlab(TeX('$\\gamma$')) + ylab('') +
  labs(color = 'Legend') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size))

# Add labels
plot5_row_form <- plot5 + shared_color_scale_row3 + shared_linetype_scale_row3 +
  guides(color = guide_legend(override.aes = list(
    linetype = c('solid', 'solid', 'solid', 'dashed'),
    linewidth = c(legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thin)
  )),
  linetype = 'none')
# plot5 <- plot5 + shared_color_scale + shared_linetype_scale +
#   guides(color = guide_legend(override.aes = list(
#     linetype = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed'),
#     linewidth = 0.6
#   )),
#   linetype = 'none')

plot5_row_form

# ------------
# -- Plot 6 --
# ------------
# Bernoulli xi
plot_df <- data.frame(eps_gamma = gamma_check,
                      avg_xi = sim6_avg_xi,
                      avg_xi_hat = sim6_avg_xi_hat,
                      avg_theta = sim6_avg_theta) %>%
  pivot_longer(c('avg_xi', 'avg_xi_hat', 'avg_theta')) %>%
  complete(name = unified_plot_levels_row3)
plot_df_ci_lower <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim6_avg_xi_lower) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row3)
plot_df_ci_upper <- data.frame(eps_gamma = gamma_check,
                               ci_bounds_xi = sim6_avg_xi_upper) %>%
  pivot_longer(c('ci_bounds_xi')) %>%
  complete(name = unified_plot_levels_row3)

# Give consistent factor levels to each of the plot_dfs
plot_df$name <- factor(plot_df$name, levels = unified_plot_levels)
plot_df_ci_lower$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row3)
plot_df_ci_upper$name <- factor(plot_df_ci_lower$name, levels = unified_plot_levels_row3)

plot6 <- ggplot() +
  geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thick_linewidth, alpha = line_alpha) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name)) +
  # geom_line(data = plot_df, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  geom_line(data = plot_df_ci_lower, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  geom_line(data = plot_df_ci_upper, aes(x = eps_gamma, y = value, color = name, linetype = name), na.rm = TRUE, linewidth = thin_linewidth, alpha = 0.9) +
  # geom_line(data = plot_df_dud, aes(x = eps_gamma, y = value, color = name), linewidth = 0.8, alpha = 0.9) +
  xlab(TeX('$\\gamma$')) + ylab('') +
  labs(color = 'Legend') +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = axis.text.x.size),
        axis.text.y = element_text(size = axis.text.y.size),
        axis.title.x = element_text(size = axis.title.x.size),
        axis.title.y = element_text(size = axis.title.y.size),
        legend.title = element_text(size = legend.title.size),
        legend.text = element_text(size = legend.text.size))

# Add labels
plot6_row_form <- plot6 + shared_color_scale_row3 + shared_linetype_scale_row3 +
  guides(color = guide_legend(override.aes = list(
    linetype = c('solid', 'solid', 'solid', 'dashed'),
    linewidth = c(legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thick, legend.linewidth.thin)
  )),
  linetype = 'none')
# plot6 <- plot6 + shared_color_scale + shared_linetype_scale +
#   guides(color = guide_legend(override.aes = list(
#     linetype = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed'),
#     linewidth = 0.6
#   )),
#   linetype = 'none')

plot6_row_form

# ===================================
# == Combine plots into one figure ==
# ===================================

# ----------------------------------------
# - Put all six plots onto a single grid -
# ----------------------------------------

final_plot <- (plot1 + plot2 + plot3 + plot4 + plot5 + plot6) +
  plot_layout(ncol = 2, guides = "collect")
final_plot

width <- 20
ggsave('figures/bernoulli_demonstration.pdf', plot = final_plot,
       device = 'pdf', width = width,
       units = 'cm')

# ----------------------------------
# - Create three rows of two plots -
# ----------------------------------
plot_row1 <- (plot1_row_form + plot2_row_form) +
  plot_layout(ncol = 2, guides = "collect")
plot_row1

plot_row2 <- (plot3_row_form + plot4_row_form) +
  plot_layout(ncol = 2, guides = "collect")
plot_row2

plot_row3 <- (plot5_row_form + plot6_row_form) +
  plot_layout(ncol = 2, guides = "collect")
plot_row3

width <- 20
height <- width * (3.75 / 10)

ggsave('figures/bernoulli_demonstration_row1.pdf', plot = plot_row1,
       device = 'pdf', width = width, height = height,
       units = 'cm')
ggsave('figures/bernoulli_demonstration_row2.pdf', plot = plot_row2,
       device = 'pdf', width = width, height = height,
       units = 'cm')
ggsave('figures/bernoulli_demonstration_row3.pdf', plot = plot_row3,
       device = 'pdf', width = width, height = height,
       units = 'cm')

