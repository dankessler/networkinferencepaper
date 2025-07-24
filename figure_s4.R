# =================================== #
# == R code to replicate Figure S4 == #
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

set.seed(14) # This seed makes the visualization look pretty good
dolphins_viz <- graph_from_adjacency_matrix(dolphins, mode = c('undirected'),
                                        weighted = FALSE, diag = FALSE)
layout <- layout_with_graphopt(dolphins_viz, charge = 0.05, niter = 500)

dolphins_viz_plot <- plot(dolphins_viz, vertex.size = 13, vertex.color = 'seashell1',
                      vertex.label.font = 2, vertex.label.color = 'black',
                      vertex.label.family = 'Helvetica',
                      vertex.label.cex = 1.2, edge.width=1.4, edge.color='black',
                      margin = 0, layout = layout)


# (Note for reader: I had a difficult time creating a line of code to save the
# plot at the exact margins and aspect ratio that I was hoping for, so I
# used a screenshot tool with the plot stretched out large in RStudio so
# that I could get a good resolution.)
