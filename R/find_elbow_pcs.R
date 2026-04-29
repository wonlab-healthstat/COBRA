#' Find the Optimal Number of PCs via the Elbow Method
#'
#' Performs truncated PCA and identifies the elbow point using the
#' Kneedle algorithm. Variance-explained values are normalised to
#' \eqn{[0, 1]} and the PC with maximum perpendicular distance from the
#' line connecting the first and last values is selected.
#' A minimum of 5 PCs is enforced.
#'
#' @param data_mat Numeric matrix (cells-by-genes or genes-by-cells).
#' @param max_pcs Maximum PCs to compute (default 50).
#' @param input_orientation \code{"cell_by_gene"} (default) or
#'   \code{"gene_by_cell"}.
#'
#' @return A list with elements \code{n_pcs}, \code{var_explained},
#'   \code{cum_var}, \code{pca_res}, and \code{elbow_plot}.
#'
#' @examples
#' \dontrun{
#' elbow <- find_elbow_pcs(expr, max_pcs = 30)
#' elbow$n_pcs
#' elbow$elbow_plot
#' }
#' @export
find_elbow_pcs <- function(data_mat,
                           max_pcs = 50,
                           input_orientation = "cell_by_gene") {

  if (input_orientation == "gene_by_cell") data_mat <- t(data_mat)

  max_pcs <- min(max_pcs, ncol(data_mat) - 1, nrow(data_mat) - 1)

  pca_res       <- irlba::prcomp_irlba(data_mat, n = max_pcs,
                                       center = TRUE, scale. = FALSE)
  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
  cum_var       <- cumsum(var_explained)

  if (length(var_explained) >= 3) {
    n      <- length(var_explained)
    x_norm <- (seq_len(n) - 1) / (n - 1)
    y_norm <- (var_explained - var_explained[n]) /
              (var_explained[1] - var_explained[n])
    distances <- abs(x_norm + y_norm - 1) / sqrt(2)
    elbow_pc  <- max(which.max(distances), 5L)
  } else {
    elbow_pc <- max_pcs
  }

  elbow_df <- data.frame(PC = seq_len(max_pcs),
                         Var_Explained = var_explained,
                         Cum_Var = cum_var)

  p_elbow <- ggplot2::ggplot(elbow_df,
                             ggplot2::aes(x = PC, y = Var_Explained)) +
    ggplot2::geom_line(linewidth = 0.8, color = COBRA_COLS["blue"]) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_vline(xintercept = elbow_pc, linetype = "dashed",
                        color = COBRA_COLS["red"], linewidth = 0.7) +
    ggplot2::annotate("text", x = elbow_pc + 1.5,
                      y = max(var_explained) * 0.8,
                      label = paste0("Elbow = ", elbow_pc),
                      hjust = 0, size = 3.5,
                      color = COBRA_COLS["red"]) +
    ggplot2::labs(x = "Principal Component",
                  y = "Proportion of Variance Explained",
                  title = "Elbow Plot for PC Selection") +
    BASE_THEME

  list(n_pcs = elbow_pc, var_explained = var_explained,
       cum_var = cum_var, pca_res = pca_res, elbow_plot = p_elbow)
}
