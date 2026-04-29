# COBRA <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

**COrrection of BAtch effect** for single-cell RNA-seq data.

COBRA removes unwanted batch variation from expression matrices while
preserving biological signal, using orthogonal projection of
batch-related design matrix columns.

## Features

- **Automatic pseudo-cell estimation** — when cell-type labels are
  unavailable, COBRA infers them via iterative K-means clustering on
  PCA-reduced space with **elbow-based PC selection** (Kneedle
  algorithm).
- **Stratified correction** — handles batch-specific (isolated) cell
  types separately from shared cell types.
- **Sparse output** — returns `dgCMatrix` objects for memory efficiency.
- **Diagnostic plots** — elbow plot and convergence history are returned
  for inspection.

## Installation

```r
# Install from GitHub
devtools::install_github("yourusername/COBRA")
```

## Quick Start

```r
library(COBRA)

# With known cell types
adj <- BatchCorrect(expr_genes_by_cells, meta,
                    batch.var = "batch",
                    cell.var  = "CellType")

# Without cell types (automatic estimation)
res <- BatchCorrect(expr_genes_by_cells, meta,
                    batch.var = "batch")

res$adj            # corrected expression matrix
res$pseudocell     # estimated cell-type labels
res$diagnostics$elbow_res$elbow_plot   # elbow diagnostic plot
```

## Citation

If you use COBRA in your research, please cite:

> Seo, S. *et al.* (2026). COBRA: batch effect correction for
> single-cell RNA-seq via orthogonal projection. *In preparation*.
