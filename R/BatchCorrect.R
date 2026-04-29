#' Batch Effect Correction for Single-cell RNA-seq Data
#'
#' Removes unwanted batch variation from an expression matrix while
#' preserving biological signal via orthogonal projection.
#' When cell-type annotations are provided they are incorporated into
#' the model; otherwise pseudo-cell types are estimated automatically
#' via \code{\link{EstimatePseudoCell}}.
#'
#' @param data.mat Expression matrix. Genes-by-cells when
#'   \code{reduce = FALSE}; cells-by-genes when \code{reduce = TRUE}.
#' @param meta.mat Data frame with cell metadata containing at least
#'   the batch column.
#' @param batch.var Column name for batch (default \code{"batch"}).
#' @param cell.var Column name for cell type, or \code{NULL} for
#'   automatic estimation (default \code{NULL}).
#' @param coi.var Additional covariates to preserve (default NULL).
#' @param stratify Stratify correction by cell type (default FALSE).
#' @param reduce TRUE if \code{data.mat} is already reduced
#'   (default FALSE).
#' @param sparse.tol Sparsification threshold (default 1e-06).
#' @param max_pcs_elbow Max PCs for elbow selection (default 50).
#'
#' @return When \code{cell.var} is provided, a sparse corrected matrix.
#'   When \code{cell.var = NULL}, a list with \code{adj} (corrected
#'   matrix), \code{pseudocell} (labels), and \code{diagnostics}.
#'
#' @examples
#' \dontrun{
#' # With known cell types
#' adj <- BatchCorrect(expr, meta, batch.var = "batch",
#'                     cell.var = "CellType")
#'
#' # Without cell types (automatic)
#' res <- BatchCorrect(expr, meta, batch.var = "batch")
#' res$adj            # corrected matrix
#' res$pseudocell     # estimated labels
#' }
#'
#' @seealso \code{\link{EstimatePseudoCell}}, \code{\link{find_elbow_pcs}}
#' @export
BatchCorrect <- function(data.mat,
                         meta.mat,
                         batch.var     = "batch",
                         cell.var      = NULL,
                         coi.var       = NULL,
                         stratify      = FALSE,
                         reduce        = FALSE,
                         sparse.tol    = 1e-06,
                         max_pcs_elbow = 50) {

  dat <- if (reduce) as.matrix(data.mat)
         else Rfast::transpose(as.matrix(data.mat))

  # Pseudo-cell estimation when cell labels are unknown
  pseudo_diag <- NULL

  if (is.null(cell.var)) {
    pca_input  <- if (reduce) data.mat else t(data.mat)
    pseudo_res <- EstimatePseudoCell(
      pca_input, batch = factor(meta.mat[, batch.var]),
      max_pcs_elbow = max_pcs_elbow
    )
    meta.mat$cell <- pseudo_res$cluster
    orig.cell.var <- NULL
    cell.var      <- "cell"
    pseudo_diag   <- pseudo_res
    cat(paste0("Number of clusters: ",
               length(unique(meta.mat$cell)), "\n"))
    cat(paste0("Number of PCs used: ", pseudo_res$n_pcs, "\n"))
  } else {
    orig.cell.var <- cell.var
  }

  # Batch effect removal
  if (!stratify) {
    meta.mat[, batch.var] <- factor(meta.mat[, batch.var])

    if (length(unique(meta.mat[, cell.var])) == 1) {
      form <- if (is.null(coi.var)) paste0("~", batch.var)
              else paste0("~", batch.var, "*(",
                          paste0(coi.var, collapse = "+"), ")")
    } else {
      meta.mat[, cell.var] <- factor(meta.mat[, cell.var])
      form <- if (is.null(coi.var)) paste0("~", batch.var, "*(", cell.var, ")")
              else paste0("~", batch.var, "*", cell.var, "*(",
                          paste0(coi.var, collapse = "+"), ")")
    }

    design    <- sparse.model.matrix(as.formula(form), data = meta.mat)
    batch_loc <- grep(batch.var, colnames(design))
    adj       <- adj_ortho(dat, design, batch_loc, reduce)

  } else {
    adj <- matrix(NA_real_, nrow(dat), ncol(dat))

    cell_prop     <- table(meta.mat[, cell.var],
                           meta.mat[, batch.var]) /
                     c(table(meta.mat[, cell.var]))
    cell_prop_idx <- apply(cell_prop, 1, function(x) sum(x > 0.99))
    cell_even     <- names(which(cell_prop_idx == 0))

    sub_cell_loc <- which(meta.mat[, cell.var] %in% cell_even)
    sub_meta     <- meta.mat[sub_cell_loc, ]
    sub_dat      <- Rfast::transpose(dat[sub_cell_loc, ])

    adj_wo_str <- BatchCorrect(data.mat = sub_dat, meta.mat = sub_meta,
                               batch.var = batch.var, cell.var = cell.var,
                               stratify = FALSE)
    if (is.list(adj_wo_str)) adj_wo_str <- adj_wo_str[[1]]
    adj[sub_cell_loc, ] <- as.matrix(adj_wo_str)

    if (length(sub_cell_loc) != nrow(adj)) {
      form <- if (is.null(coi.var)) paste0("~", batch.var, "+", cell.var)
              else paste0("~", batch.var, "+", cell.var, "+", coi.var)
      design    <- sparse.model.matrix(as.formula(form), data = meta.mat)
      batch_loc <- grep(batch.var, colnames(design))
      tmp_adj   <- adj_ortho(dat, design, batch_loc, reduce)
      adj[-sub_cell_loc, ] <- tmp_adj[-sub_cell_loc, ]
    }
  }

  adj[abs(adj) < sparse.tol] <- 0
  adj <- Matrix(adj, sparse = TRUE)

  if (is.null(orig.cell.var)) {
    list(adj = adj, pseudocell = meta.mat$cell, diagnostics = pseudo_diag)
  } else {
    adj
  }
}
