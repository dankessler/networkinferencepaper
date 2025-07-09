# ==================================================== #
# == R code to replicate the left panel of Figure 8 == #
# == Author: Ethan Ancell                           == #
# ==================================================== #

library("networkinference")
library("latex2exp")
library("tidyverse")
library("igraph")

ggplot2::theme_set(theme_minimal())

# ===================================================== #
# == Manually load the Zachary's karate club dataset == #
# ===================================================== #

n <- 34
K <- 2

zachary <- matrix(0, nrow = n, ncol = n)

# List of connections (on lower-triangular portion)
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

# Fill in upper triangular portion based upon connections
for (connection in connections) {
  zachary[connection[2], connection[1]] <- 1 # Swap order so it gets filled on upper triangular
}

# Set the lower diagonal of the matrix to also be the same
zachary_upper_diagonal <- zachary
zachary <- zachary + t(zachary)

# ========== #
# == Plot == #
# ========== #

zachary_visual_df <- reshape2::melt(zachary)
colnames(zachary_visual_df) <- c('Row', 'Column', 'Value')

zachary_visual <- ggplot(zachary_visual_df) +
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

zachary_visual
ggsave('figures/zachary/zachary_cartoon_adjacency_matrix.pdf', zachary_visual, device = 'pdf',
       width = 7, height = 7)
# TODO: this saved title is not the same as it is currently listed in Overleaf so be careful!
