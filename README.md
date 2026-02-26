# Plasmid Replicon and Resistance Gene Association Analysis

This repository contains the R scripts used to generate the figures presented in the manuscript submitted to *Current Microbiology*.

## Overview

Two main analyses are included:

### 1. Correspondence Analysis (CA)

Evaluates the association between plasmid incompatibility groups and resistance categories (antibiotic, metal, and biocide resistance genes).

Output:
- Two-dimensional CA plot (Dim 1 × Dim 2)
- Variance explained by principal axes

### 2. Pairwise Gene Association Heatmap

Displays positive and negative associations between resistance genes using a signed score:

signed_score = sign(log(OR)) × -log10(q_FDR)

- Positive associations: red
- Negative associations: blue
- Intensity proportional to normalized score magnitude

---

## Required R Version

R ≥ 4.0.0

---

## Required Packages

- data.table
- ggplot2
- ca
- ggrepel (optional)

Install if necessary:

```r
install.packages(c("data.table", "ggplot2", "ca", "ggrepel"))
```

---

## Input Files

Place the following files inside a folder named `data/`:

- contingency_inc_vs_category.csv
- gene_pairwise_association_TOP.csv

---

## How to Run

Set your working directory and define:

```r
out_dir <- "data"
```

Then run:

```r
source("01_MCA_analysis.R")
source("02_Heatmap_pairwise_association.R")
```

---

## Reproducibility

Session information can be generated using:

```r
sessionInfo()
```

---

## Author

Hannay Souza
