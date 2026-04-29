#' COBRA: COrrection of BAtch Effect for Single-Cell RNA-Seq Data
#'
#' @description
#' Removes unwanted batch variation from single-cell RNA-seq expression
#' matrices while preserving biological signal. When cell-type labels are
#' unavailable, pseudo-cell types are estimated automatically via
#' iterative clustering with elbow-based PC selection.
#'
#' @import Matrix
#' @importFrom ClusterR Optimal_Clusters_KMeans KMeans_arma predict_KMeans
#' @importFrom Rfast transpose
#' @importFrom MASS ginv
#' @importFrom irlba prcomp_irlba
#' @importFrom stats as.formula lm.fit model.matrix
#' @docType package
#' @name COBRA-package
"_PACKAGE"

# Suppress R CMD check NOTEs for ggplot2 aes variables
utils::globalVariables(c("PC", "Var_Explained"))

# Internal constants
COBRA_COLS <- c(blue = "#2E86AB", red = "#E74C3C", grey = "#7F8C8D")

BASE_THEME <- ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(size = 12, face = "bold"),
    axis.title       = ggplot2::element_text(size = 10),
    axis.text        = ggplot2::element_text(size = 9),
    panel.grid.minor = ggplot2::element_blank()
  )
