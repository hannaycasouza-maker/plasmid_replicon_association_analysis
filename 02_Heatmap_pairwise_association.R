# ==========================================================
# Script: Pairwise Gene Association Heatmap
# Author: Hannay Souza
# Journal: Current Microbiology
# Description:
# Generates heatmap of pairwise gene associations
# based on signed log(OR) × -log10(q_FDR)
# ==========================================================

rm(list = ls())

library(data.table)
library(ggplot2)

out_dir <- "data"

assoc_path <- file.path(out_dir, "gene_pairwise_association_TOP.csv")
stopifnot(file.exists(assoc_path))

assoc <- data.table::fread(assoc_path)

if (!("i" %in% names(assoc) && "j" %in% names(assoc))) {
  if (all(c("gene_i","gene_j") %in% names(assoc))) {
    data.table::setnames(assoc, c("gene_i","gene_j"), c("i","j"))
  } else {
    stop("Gene pair columns not recognized.")
  }
}

if (!("signed_score" %in% names(assoc))) {
  stop("Column 'signed_score' not found.")
}

labs <- sort(unique(c(assoc$i, assoc$j)))

M <- matrix(0,
            nrow = length(labs),
            ncol = length(labs),
            dimnames = list(labs, labs))

for (k in seq_len(nrow(assoc))) {
  M[assoc$i[k], assoc$j[k]] <- assoc$signed_score[k]
  M[assoc$j[k], assoc$i[k]] <- assoc$signed_score[k]
}

hm <- as.data.frame(as.table(M))
names(hm) <- c("gene_x","gene_y","score")

max_abs <- max(abs(hm$score), na.rm = TRUE)

p_heat <- ggplot(hm, aes(gene_x, gene_y, fill = score)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "sign(log(OR)) × -log10(q_FDR)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Pairwise Association Among Resistance Genes",
    x = NULL,
    y = NULL
  )

print(p_heat)
