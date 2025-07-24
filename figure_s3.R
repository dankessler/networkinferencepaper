# =================================== #
# == R code to replicate Figure S3 == #
# == Author: Ethan Ancell          == #
# =================================== #

library("manynet") # Where dolphins dataset lives
library("networkinference")
library("latex2exp")
library("tidyverse")
library("igraph")

ggplot2::theme_set(theme_minimal())

# ====================================
# == Load dolphins adjacency matrix ==
# ====================================

data("ison_dolphins")
dolphins_igraph <- as_igraph(ison_dolphins)

# Get adjacency matrix
dolphins <- igraph::as_adjacency_matrix(dolphins_igraph, sparse = FALSE)

# Remove names because they crowd the screen
colnames(dolphins) <- 1:NCOL(dolphins)
rownames(dolphins) <- 1:NROW(dolphins)

# ========== #
# == Plot == #
# ========== #

dolphins_visual_df <- reshape2::melt(dolphins)
colnames(dolphins_visual_df) <- c('Row', 'Column', 'Value')

dolphins_visual <- ggplot(dolphins_visual_df) +
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

dolphins_visual
ggsave('figures/dolphins_cartoon_adjacency_matrix.pdf', dolphins_visual, device = 'pdf',
       width = 7, height = 7)
