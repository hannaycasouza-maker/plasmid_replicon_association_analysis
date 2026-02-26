# ==========================================================
# Script: Correspondence Analysis (CA)
# Author: Hannay Souza
# Journal: Current Microbiology
# Description:
# Performs correspondence analysis between plasmid
# incompatibility groups and resistance categories.
# ==========================================================

rm(list = ls())

library(data.table)
library(ggplot2)
library(ca)

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# Define data directory
out_dir <- "data"

cont_path <- file.path(out_dir, "contingency_inc_vs_category.csv")
stopifnot(file.exists(cont_path))

Cont <- data.table::fread(cont_path)

label_col <- if ("incompatibility_group" %in% names(Cont)) {
  "incompatibility_group"
} else {
  cand <- names(Cont)[vapply(Cont, function(x) is.character(x) || is.factor(x), logical(1))]
  if (length(cand)) cand[1] else names(Cont)[1]
}

rn <- as.character(Cont[[label_col]])
Cont[[label_col]] <- NULL

for (nm in names(Cont)) Cont[[nm]] <- suppressWarnings(as.numeric(Cont[[nm]]))
Cont[is.na(Cont)] <- 0

mat <- as.matrix(Cont)
rownames(mat) <- rn

keep_r <- rowSums(mat) > 0
keep_c <- colSums(mat) > 0
mat <- mat[keep_r, keep_c, drop = FALSE]

if (nrow(mat) < 2 || ncol(mat) < 2) stop("Matrix insufficient for CA.")

fit <- ca::ca(mat)

rows <- data.frame(
  Dim1 = fit$rowcoord[,1],
  Dim2 = fit$rowcoord[,2],
  label = rownames(fit$rowcoord)
)

cols <- data.frame(
  Dim1 = fit$colcoord[,1],
  Dim2 = fit$colcoord[,2],
  label = rownames(fit$colcoord)
)

inertia <- (fit$sv^2) / sum(fit$sv^2)

ttl <- sprintf(
  "Correspondence Analysis (Inc group × Resistance category)\nDim1: %.1f%%  Dim2: %.1f%%",
  100 * inertia[1],
  100 * inertia[2]
)

p_ca <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey35") +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey35") +
  geom_point(data = rows, aes(Dim1, Dim2), shape = 4, size = 2.8, colour = "blue") +
  geom_point(data = cols, aes(Dim1, Dim2), shape = 22, size = 3.5,
             fill = "red", colour = "red") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  labs(title = ttl, x = "Dim 1", y = "Dim 2")

print(p_ca)
