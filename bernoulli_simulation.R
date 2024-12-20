# ============================================
# == Simulations for network thinning paper ==
# == (Binary edges)                         ==
# == Author: Ethan Ancell                   ==
# ============================================

# Libraries
library(ggplot2)
library(tidyverse)
library(nett)
library(mclust)
library(latex2exp)

# =========================
# == Simulation settings ==
# =========================

set.seed(1)

# Most comprehensive settings
# (Currently commented out)
# ---------------------------

# n_check <- c(50, 100, 200)                   # Total network size
# n_per_community_check <- c(1, 2, 10)         # Amount of nodes in each true community
# K_check <- c(2, 5, 10)                       # Guessed number of communities
# num_sim <- 500                               # Number of replications
# signal_regimes <- list(c(11, 10), c(15, 10)) # Signal regimes for mean matrix
# eps_check <- c(0.2, 0.5, 0.7)


# Feel free to adjust any of these
# settings on the fly for whatever you're doing
# ---------------------------------------------

# -------------
# -- Table 1 --
# -------------

set.seed(1)
num_sim <- 5000 # Number of replications

n_check <- c(100, 200, 500) # Total network size
K_check <- c(2, 5, 10) # Guessed number of communities
K_true_check <- c(5)
eps_check <- c(0.25)
signal_regimes <- list(c(0.75, 0.5)) # Signal regimes for mean matrix

alpha <- 0.10
use_random_u <- FALSE

# --------------
# -- Figure 1 --
# --------------

set.seed(1)
num_sim <- 5000 # Number of replications

n_check <- c(200) # Total network size
K_check <- c(2, 5, 10) # Guessed number of communities
K_true_check <- c(5)
eps_check <- c(0.25)
signal_regimes <- list(c(0.75, 0.5)) # Signal regimes for mean matrix

alpha <- 0.10
use_random_u <- FALSE

# --------------
# -- Figure 2 --
# --------------

set.seed(1)
num_sim <- 5000 # Number of replications

n_check <- c(200) # Total network size
K_check <- c(5) # Guessed number of communities
K_true_check <- c(5)
eps_check <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)
signal_regimes <- list(c(0.75, 0.55), c(0.75, 0.50), c(0.75, 0.45), c(0.75, 0.40), c(0.75, 0.35)) # Signal regimes for mean matrix

alpha <- 0.10
use_random_u <- FALSE

# =======================
# == Set up simulation ==
# =======================

# Useful for bookkeeping and plotting: which dimensions correspond to which
# things that we are looping through.
looping_dimensions <- c('n', 'K_true', 'K', 'Signal Regime')

# Unpack rho_1 and rho_2 from the signal regimes
rho_1 <- signal_regimes[[1]][1] # Assumes that this stays the same across all of them
rho_2_check <- rep(NA, length(signal_regimes))
for (i in 1:length(signal_regimes)) {
  rho_2_check[i] <- signal_regimes[[i]][2]
}
rho_1_minus_rho_2 <- rho_1 - rho_2_check

# Where results are stored
# ------------------------

# Estimates, targets, variances, etc.
# (saving all of this makes it easier to calculate some things without
# re-running the entire experiment)
estimates_thinning <- array(0, dim = c(length(n_check), 
                                       length(K_true_check),
                                       length(K_check),
                                       length(signal_regimes),
                                       length(eps_check),
                                       num_sim))
targets_thinning <- array(0, dim = c(length(n_check), 
                                     length(K_true_check),
                                     length(K_check),
                                     length(signal_regimes),
                                     length(eps_check),
                                     num_sim))
variances_thinning <- array(0, dim = c(length(n_check), 
                                       length(K_true_check),
                                       length(K_check),
                                       length(signal_regimes),
                                       length(eps_check),
                                       num_sim))

estimates_naive <- array(0, dim = c(length(n_check), 
                                    length(K_true_check),
                                    length(K_check),
                                    length(signal_regimes),
                                    length(eps_check),
                                    num_sim))
targets_naive <- array(0, dim = c(length(n_check), 
                                  length(K_true_check),
                                  length(K_check),
                                  length(signal_regimes),
                                  length(eps_check),
                                  num_sim))
variances_naive <- array(0, dim = c(length(n_check), 
                                    length(K_true_check),
                                    length(K_check),
                                    length(signal_regimes),
                                    length(eps_check),
                                    num_sim))
rand_results <- array(0, dim = c(length(n_check), 
                                 length(K_true_check),
                                 length(K_check),
                                 length(signal_regimes),
                                 length(eps_check),
                                 num_sim))

# Coverage results of (1-alpha) confidence intervals
covered_thinning <- array(0, dim = c(length(n_check), 
                                     length(K_true_check),
                                     length(K_check),
                                     length(signal_regimes),
                                     length(eps_check)))
covered_naive <- array(0, dim = c(length(n_check), 
                                  length(K_true_check),
                                  length(K_check),
                                  length(signal_regimes),
                                  length(eps_check)))

# ==============
# == Simulate ==
# ==============

# Sample size
for (n_index in 1:length(n_check)) {
  n <- n_check[n_index]
  cat(paste0('(n = ', n, ')\n'))
  
  # True communities
  for (K_true_index in 1:length(K_true_check)) {
    # n_per_community <- n_per_community_check[n_per_community_index]
    K_true <- K_true_check[K_true_index]
    n_per_community <- n / K_true
    cat(paste0('\t(K_true = ', K_true, ')\n'))
    
    # Number of communities used in clustering
    for (K_index in 1:length(K_check)) {
      K <- K_check[K_index]
      cat(paste0('\t\t(K = ', K, ')\n'))
      
      # Signal regimes in mean matrix
      for (signal_regime_index in 1:length(signal_regimes)) {
        rho1 <- signal_regimes[[signal_regime_index]][1]
        rho2 <- signal_regimes[[signal_regime_index]][2]
        cat(paste0('\t\t\t(rho1 = ', rho1, ', rho2 = ', rho2, ')\n'))
        
        # Different values of epsilon in thinning
        for (eps_index in 1:length(eps_check)) {
          eps <- eps_check[eps_index]
          cat(paste0('\t\t\t\t(eps = ', eps, ')\n'))
          
          # -------------------------------------------------
          # -- Actual simulation part (non just for-loops) --
          # -------------------------------------------------
          
          # Create mean matrices
          Theta_true <- matrix(rep(rho2, K_true^2), nrow = K_true)
          diag(Theta_true) <- rho1
          M <- matrix(rep(NA, n^2), nrow = n)
          for (k in 1:K_true) {
            for (l in 1:K_true) {
              M[((k-1)*n_per_community + 1):(k*n_per_community),
                ((l-1)*n_per_community + 1):(l*n_per_community)] <- Theta_true[k, l]
            }
          }
          # "True" community membership
          z_true <- rep(1:K_true, each = n_per_community)
          
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
            A <- matrix(rbinom(n = n^2, size = 1, prob = as.vector(M)), nrow = n)
            # A_2 <- matrix(rbinom(n = n^2, size = 1, prob = as.vector(M)), nrow = n)
            
            # Fission A
            Zfission <- matrix(rbinom(n = n^2, size = 1, prob = eps), nrow = n)
            A_tr <- A * (1 - Zfission) + (1 - A) * Zfission
            
            # Cluster (thinning)
            z_hat_initial <- nett::spec_clust(A_tr, K = K)
            z_hat <- matrix(rep(NA, n*K), nrow = n)
            for (i in 1:K) {
              z_hat[, i] <- 1 * (z_hat_initial == i)
            }
            n_hat <- apply(z_hat, 2, sum)
            
            # Check RAND index agreement between true clusters and non true clusters
            rand_results[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <- 
              mclust::adjustedRandIndex(z_hat_initial, z_true)
            
            # Cluster (naive)
            # z_hat_initial_naive <- nett::spec_clust(A_2, K = K)
            z_hat_initial_naive <- nett::spec_clust(A, K = K)
            z_hat_naive <- matrix(rep(NA, n*K), nrow = n)
            for (i in 1:K) {
              z_hat_naive[, i] <- 1 * (z_hat_initial_naive == i)
            }
            n_hat_naive <- apply(z_hat_naive, 2, sum)
            
            # ------------------
            # -- EXPERIMENTAL --
            # ------------------
            
            # Cmask <- (eps / (1-eps))^(2*A_tr - 1) # The constant epsilon/(1-epsilon) masking thing
            # Tmat <- M / (M + (1-M)*Cmask)
            # 
            # # Draw a faux copy of A conditional on Atr
            # A_faux <- matrix(rbinom(n = n^2, size = 1, prob = as.vector(Tmat)), nrow = n)
            # A <- A_faux
            
            # ----------------------------------
            # -- Back to the non-experimental --
            # ----------------------------------
            
            # Target and estimator (thinning)
            Cmask <- (eps / (1-eps))^(2*A_tr - 1) # The constant epsilon/(1-epsilon) masking thing
            Tmat <- M / (M + (1-M)*Cmask) # Calculate T as a function of M
            NN_inv <- diag(1 / diag(t(z_hat) %*% z_hat))
            estimate_blocky <- NN_inv %*% t(z_hat) %*% A %*% z_hat %*% NN_inv
            estimate <- t(u) %*% as.vector(NN_inv %*% t(z_hat) %*% A %*% z_hat %*% NN_inv)
            target <- t(u) %*% as.vector(NN_inv %*% t(z_hat) %*% Tmat %*% z_hat %*% NN_inv)
            
            # Target and estimator (naive)
            NN_inv_naive <- diag(1 / diag(t(z_hat_naive) %*% z_hat_naive))
            estimate_blocky_naive <- NN_inv_naive %*% t(z_hat_naive) %*% A %*% z_hat_naive %*% NN_inv_naive
            estimate_naive <- t(u) %*% as.vector(NN_inv_naive %*% t(z_hat_naive) %*% A %*% z_hat_naive %*% NN_inv_naive)
            target_naive <- t(u) %*% as.vector(NN_inv_naive %*% t(z_hat_naive) %*% M %*% z_hat_naive %*% NN_inv_naive)
            
            # Estimated variance of estimator (thinning and naive)
            Delta_n <- matrix(rep(NA, K^2), nrow = K)
            Delta_n_naive <- matrix(rep(NA, K^2), nrow = K)
            for (k in 1:K) {
              for (l in 1:K) {
                # Thinning
                k_nodes <- which(z_hat_initial == k)
                l_nodes <- which(z_hat_initial == l)
                
                # Jn <- sum((A[k_nodes, l_nodes] / (estimate_blocky[k, l]^2)) + 
                #             ((1-A[k_nodes, l_nodes]) / ((1 - estimate_blocky[k, l])^2)))
                # Jn_true <- sum((Tmat[k_nodes, l_nodes] / (estimate_blocky[k, l]^2)) +
                #                  ((1 - Tmat[k_nodes, l_nodes]) / ((1 - estimate_blocky[k, l])^2)))
                # Vn_true <- sum(Tmat[k_nodes, l_nodes]*(1 - Tmat[k_nodes, l_nodes])) / (estimate_blocky[k, l] * (1 - estimate_blocky[k, l]))^2
                
                # TODO:
                # Jn <- sum(Tmat[k_nodes, l_nodes] * (1 - Tmat[k_nodes, l_nodes]))
                # Delta_n[k, l] <- sum(Tmat[k_nodes, l_nodes] * (1 - Tmat[k_nodes, l_nodes])) / (length(k_nodes) * length(l_nodes))^2
                Delta_n[k, l] <- (estimate_blocky[k, l] * (1 - estimate_blocky[k, l])) / (length(k_nodes) * length(l_nodes))
                
                #Delta_n[k, l] <- 1 / Jn
                # Delta_n[k, l] <- Vn_true / Jn^2
                
                # Naive
                k_nodes_naive <- which(z_hat_initial_naive == k)
                l_nodes_naive <- which(z_hat_initial_naive == l)
                
                # Jn_naive <- sum((A[k_nodes_naive, l_nodes_naive] / estimate_blocky_naive[k, l]^2) + 
                #                   ((1-A[k_nodes_naive, l_nodes_naive]) / (1 - estimate_blocky_naive[k, l])^2))
                
                Delta_n_naive[k, l] <- (estimate_blocky_naive[k, l] * (1 - estimate_blocky_naive[k, l])) / (length(k_nodes_naive) * length(l_nodes_naive))
              }
            }
            estimate_var <- t(u) %*% diag(as.vector(Delta_n)) %*% u
            estimate_var_naive <- t(u) %*% diag(as.vector(Delta_n_naive)) %*% u
            
            # Construct a confidence interval for the target of inference
            margin_of_error <- qnorm(1 - alpha / 2) * sqrt(estimate_var)
            margin_of_error_naive <- qnorm(1 - alpha / 2) * sqrt(estimate_var_naive)
            
            # Check coverage (thinning and naive)
            covered_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index] + 
              as.integer((estimate - margin_of_error <= target) & (target <= estimate + margin_of_error))
            
            covered_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index] <-
              covered_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index] + 
              as.integer((estimate_naive - margin_of_error_naive <= target_naive) & (target_naive <= estimate_naive + margin_of_error_naive))
            
            # Save everything else
            # --------------------
            estimates_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate
            targets_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target
            variances_thinning[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var
            estimates_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_naive
            targets_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              target_naive
            variances_naive[n_index, K_true_index, K_index, signal_regime_index, eps_index, rep] <-
              estimate_var_naive
          }
        }
      }
    }
  }
}

covered_thinning <- covered_thinning / num_sim
covered_naive <- covered_naive / num_sim

covered_thinning
covered_naive

# ==================================
# == Plotting and viewing results ==
# ==================================

# Table 1 - Checking coverage of thinning and naive
# -------------------------------------------------
plot_table_1 <- function(data) {
  
  signal_index <- 1
  K_true_index <- 1
  eps_index <- 1
  
  # Print header column
  cat("|         |")
  for (K_index in 1:length(K_check)) {
    K <- K_check[K_index]
    cat(sprintf(" K = %03d |", K))
  }
  cat('\n')
  
  # Print rows of table below header
  for (n_index in 1:length(n_check)) {
    n <- n_check[n_index]
    
    cat(sprintf("| n = %03d |", n))
    
    for (K_index in 1:length(K_check)) {
      K <- K_check[K_index]
      
      cat(sprintf("   %.2f |", 100 * data[n_index, K_true_index, K_index, signal_index, eps_index]))
    }
    cat('\n')
  }
}

cat('Naive\n')
plot_table_1(covered_naive)

cat('Thinning\n')
plot_table_1(covered_thinning)

# Figure 1 - Showing a comparison of nominal vs empirical coverage
# ----------------------------------------------------------------

checked_coverages <- seq(0, 1.0, length.out = 100)
n_coverages <- length(checked_coverages)

n_index <- 1
K_true_index <- 1
K_true <- K_true_check[K_true_index]
signal_index <- 1

# Results
figure1_coverage_thinning <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))
figure1_coverage_naive <- array(numeric(), c(length(K_check), length(eps_check), length(checked_coverages)))

for (K_index in 1:length(K_check)) {
  
  for (eps_index in 1:length(eps_check)) {
    
    for (coverage_index in 1:length(checked_coverages)) {
      
      coverage <- checked_coverages[coverage_index]
      z_quantile <- qnorm(1 - (1 - coverage) / 2)
      
      # Extract relevant slices from the results array
      targets_thinning_slice <- targets_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ]
      targets_naive_slice <- targets_naive[n_index, K_true_index, K_index, signal_index, eps_index, ]
      
      upper_bounds_thinning <- 
        estimates_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ] +
        z_quantile * sqrt(variances_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ])
      lower_bounds_thinning <- 
        estimates_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ] -
        z_quantile * sqrt(variances_thinning[n_index, K_true_index, K_index, signal_index, eps_index, ])
      
      upper_bounds_naive <- 
        estimates_naive[n_index, K_true_index, K_index, signal_index, eps_index, ] +
        z_quantile * sqrt(variances_naive[n_index, K_true_index, K_index, signal_index, eps_index, ])
      lower_bounds_naive <- 
        estimates_naive[n_index, K_true_index, K_index, signal_index, eps_index, ] -
        z_quantile * sqrt(variances_naive[n_index, K_true_index, K_index, signal_index, eps_index, ])
      
      # Store coverage
      figure1_coverage_thinning[K_index, eps_index, coverage_index] <-
        mean((lower_bounds_thinning <= targets_thinning_slice) & (targets_thinning_slice <= upper_bounds_thinning))
      figure1_coverage_naive[K_index, eps_index, coverage_index] <-
        mean((lower_bounds_naive <= targets_naive_slice) & (targets_naive_slice <= upper_bounds_naive))
      
    }
    
  }
  
}

# Stack all this information in a data frame that can be used for plotting
plot_df <- data.frame()
for (K_index in 1:length(K_check)) {
  for (eps_index in 1:length(eps_check)) {
    
    # Create a method name for thinning / naive that includes K = ---
    # thinning_method_name <- paste0('Thinning (K = ', K_check[K_index], ')')
    # naive_method_name <- paste0('Naive (K = ', K_check[K_index], ')')
    
    # Thinning
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_thinning[K_index, eps_index, ],
                                method = rep('Fission', n_coverages),
                                # method_print = rep(thinning_method_name, n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages)))
    
    # Naive
    plot_df <- rbind(plot_df,
                     data.frame(index = 1:n_coverages,
                                nominal_coverage = checked_coverages,
                                coverage = figure1_coverage_naive[K_index, eps_index, ],
                                method = rep('Naive', n_coverages),
                                # method_print = rep(naive_method_name, n_coverages),
                                K = rep(K_check[K_index], n_coverages),
                                Kch = as.character(rep(K_check[K_index], n_coverages)),
                                eps = rep(eps_check[eps_index], n_coverages)))
  }
}

# Create plot
legend_colors <- c('2' = 'cadetblue3', '5' = 'royalblue3', '10' = 'blue4')
figure1 <- ggplot()
for (K_index in 1:length(K_check)) {
  K <- K_check[K_index]
  
  # Subset data
  df_subset_thin <- plot_df %>% filter(method == 'Fission', K == K)
  df_subset_naive <- plot_df %>% filter(method == 'Naive', K == K)
  
  # Add it onto the plot
  figure1 <- figure1 +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Fission'), size = 0.9, data = df_subset_thin)
  figure1 <- figure1 +
    geom_line(aes(x = nominal_coverage, y = coverage, color = Kch, linetype = 'Naive'), size = 0.7, data = df_subset_naive)
}
figure1 <- figure1 +
  scale_color_manual(breaks = c('2', '5', '10'), 
                     values = legend_colors) +
  scale_linetype_manual(breaks = c('Fission', 'Naive'),
                        values = c('Fission' = 'solid', 'Naive' = 
                                     'dashed')) +
  xlab('Nominal Coverage') + ylab('Empirical Coverage') +
  coord_fixed() +
  labs(color = 'K', linetype = "Method") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

figure1
ggsave('final_simulations/figures/nominal_empirical_coverage_binary.pdf',
       plot = figure1, device = 'pdf', width = 4.5, height = 3.5)

# Figure 2 - Tradeoff of signal regime, epsilon, RAND index, CI width
# -------------------------------------------------------------------

n_index <- 1
K_true_index <- 1
K_index <- 1
alpha <- 0.10

# Figure 2a - CI width as a function of signal, separate by epsilon
# -----------------------------------------------------------------
ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_true_index, K_index, , , ])
# ci_width <- (2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_true_index, K_index, , , ])) / sqrt(estimates_thinning[n_index, K_true_index, K_index, , , ])
ci_width <- apply(ci_width, c(1, 2), mean)
eps_index_plotting <- matrix(rep(eps_check, length(signal_regimes)), ncol = length(eps_check), byrow = TRUE)
rho_2_plotting <- matrix(rep(rho_2_check, length(eps_check)), ncol = length(eps_check))
rho_1_minus_rho_2_plotting <- matrix(rep(rho_1_minus_rho_2, length(eps_check)), ncol = length(eps_check))
signal_index_plotting <- matrix(rep(1:length(signal_regimes), length(eps_check)), ncol = length(eps_check))

plot_df <- data.frame(signal_index = as.vector(signal_index_plotting),
                      rho_2 = as.vector(rho_2_plotting),
                      rho_1_minus_rho_2 = as.character(as.vector(rho_1_minus_rho_2_plotting)),
                      #signal_name = paste0('rho_2 = ', as.vector(rho_2_plotting)),
                      ci_width = as.vector(ci_width), 
                      eps = as.vector(eps_index_plotting))

legend_colors <- c('0.2' = 'gold3', '0.25' = 'darkseagreen4',
                   '0.3' = 'darkslategray4', '0.35' = 'deeppink4',
                   '0.4' = 'firebrick3')
# legend_colors = c('rho_2 = 20' = 'springgreen4', 'rho_2 = 25' = 'violetred4', 'rho_2 = 30' = 'slateblue4')
figure2a <- ggplot()
for (ri in 1:length(rho_2_check)) {
  
  # signal_name = paste0('rho2 = ', rho_2_check[ri])
  
  df_subset <- plot_df[plot_df$rho_2 == rho_2_check[ri], ]
  figure2a <- figure2a + 
    geom_line(aes(x = eps, y = ci_width, color = rho_1_minus_rho_2), size = 1.0, alpha = 0.7, data = df_subset)
}
figure2a <- figure2a +
  xlab(TeX('$\\epsilon$')) + ylab('Average CI width') +
  scale_color_manual(breaks = c('0.2', '0.25', '0.3', '0.35', '0.4'), values = legend_colors) +
  # ggtitle(TeX('Average confidence interval width as a function of $\\epsilon$')) +
  labs(color = TeX('$\\rho_1 - \\rho_2$')) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
figure2a
ggsave('final_simulations/figures/conf_width_epsilon_binary.pdf',
       plot = figure2a, device = 'pdf', width = 4.5, height = 3.5)

# Figure 2b - Tradeoff of signal regime, RAND index for clustering, separated by epsilon
# --------------------------------------------------------------------------------------
#ci_width <- 2 * qnorm(1 - alpha / 2) * sqrt(variances_thinning[n_index, K_true_index, K_index, , , ])
# ci_width <- apply(ci_width, c(1, 2), mean)
rand_slice <- apply(rand_results[n_index, K_true_index, K_index, , , ], c(1, 2), mean)
eps_index_plotting <- matrix(rep(eps_check, length(signal_regimes)), ncol = length(eps_check), byrow = TRUE)
rho_2_plotting <- matrix(rep(rho_2_check, length(eps_check)), ncol = length(eps_check))
rho_1_minus_rho_2_plotting <- matrix(rep(rho_1_minus_rho_2, length(eps_check)), ncol = length(eps_check))
signal_index_plotting <- matrix(rep(1:length(signal_regimes), length(eps_check)), ncol = length(eps_check))

legend_colors <- c('0.2' = 'gold3', '0.25' = 'darkseagreen4',
                   '0.3' = 'darkslategray4', '0.35' = 'deeppink4',
                   '0.4' = 'firebrick3')

plot_df <- data.frame(signal_index = as.vector(signal_index_plotting),
                      rho_2 = as.character(as.vector(rho_2_plotting)),
                      rho_1_minus_rho_2 = as.character(as.vector(rho_1_minus_rho_2_plotting)),
                      rand_index = as.vector(rand_slice), 
                      eps = as.vector(eps_index_plotting))
figure2b <- ggplot(plot_df) +
  geom_line(aes(x = eps, y = rand_index, color = rho_1_minus_rho_2), size = 0.9) +
  scale_color_manual(breaks = c('100', '101', '102', '103', '104'), values = legend_colors) + # breaks are wrong because I don't want a legend to appear...
  xlab(TeX('$\\epsilon')) + ylab('Adjusted RAND Index') +
  coord_fixed() +
  ylim(c(0, 1.0)) +
  #ggtitle(TeX('Classification accuracy (Adjusted RAND index) as a function of $\\epsilon$')) +
  labs(color = TeX('$\\rho_1 - \\rho_2$')) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
figure2b
ggsave('final_simulations/figures/randindex_epsilon_binary.pdf',
       plot = figure2b, device = 'pdf', width = 4.5, height = 3.5)
