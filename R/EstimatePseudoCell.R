#' Estimate Pseudo-cell Types with Elbow-based PC Selection
#'
#' When cell-type labels are unavailable, infers pseudo-cell-type
#' clusters in three stages: (1) PCA with automatic elbow-based
#' dimensionality selection, (2) K-means initial clustering with
#' silhouette-based K selection, (3) iterative refinement via
#' orthogonal batch removal until BIC converges.
#'
#' @param data_mat Numeric matrix, cells-by-genes.
#' @param batch Factor of batch labels.
#' @param ncl_min Minimum K to evaluate (default 2).
#' @param ncl_max Maximum K to evaluate (default 10).
#' @param seed Random seed (default 1).
#' @param r_random_cell Proportion of cells for initial clustering
#'   (default 0.01).
#' @param n_random_cell Minimum cells for initial clustering
#'   (default 1000).
#' @param pheno Optional phenotype covariate (default NULL).
#' @param max_pcs_elbow Maximum PCs for elbow search (default 50).
#'
#' @return A list with \code{cluster}, \code{n_pcs}, \code{history},
#'   and \code{elbow_res}.
#'
#' @examples
#' \dontrun{
#' pseudo <- EstimatePseudoCell(expr_mat, batch = factor(meta$batch))
#' table(pseudo$cluster)
#' }
#' @export
EstimatePseudoCell <- function(data_mat, batch,
                               ncl_min = 2, ncl_max = 10,
                               seed = 1,
                               r_random_cell = 0.01,
                               n_random_cell = 1000,
                               pheno = NULL,
                               max_pcs_elbow = 50) {

  set.seed(seed)
  n_sub    <- min(max(nrow(data_mat) * r_random_cell, n_random_cell),
                  nrow(data_mat))
  cell_loc <- sample(seq_len(nrow(data_mat)), n_sub)

  cat("Estimate Pseudo-cells (elbow-based PC selection)\n")

  # Step 1: Elbow-based PC selection
  elbow_res <- find_elbow_pcs(data_mat, max_pcs = max_pcs_elbow,
                              input_orientation = "cell_by_gene")
  n_pcs   <- elbow_res$n_pcs
  dat_pcs <- elbow_res$pca_res$x[, seq_len(n_pcs), drop = FALSE]
  cat(paste0("  Selected ", n_pcs, " PCs via elbow method\n"))

  # Step 2: Initial clustering
  opt_km  <- Optimal_Clusters_KMeans(dat_pcs[cell_loc, , drop = FALSE],
                                     max_clusters = ncl_min:ncl_max,
                                     plot_clusters = FALSE, verbose = FALSE,
                                     criterion = "silhouette")
  ncl     <- max(which(opt_km == max(opt_km))) + ncl_min - 1
  km      <- KMeans_arma(dat_pcs, ncl, n_iter = 100,
                         "random_subset", verbose = FALSE)
  cluster <- as.factor(predict_KMeans(dat_pcs, km))
  bic_km  <- Optimal_Clusters_KMeans(dat_pcs, max_clusters = ncl_min:ncl_max,
                                     plot_clusters = FALSE, verbose = FALSE,
                                     criterion = "BIC")

  # Step 3: Iterative refinement
  convergence_history <- data.frame(iter = integer(), k = integer(),
                                    bic = numeric(), silhouette = numeric(),
                                    accepted = logical())
  iter    <- 1
  adj_pcs <- dat_pcs

  while (iter < 20) {
    batch_prop   <- table(cluster, batch) / c(table(cluster))
    check_uneven <- ifelse(batch_prop > 0.9, 1, 0)

    if (sum(check_uneven == 1) >
        length(check_uneven) / length(unique(batch)) * 0.5) {
      design <- if (is.null(unlist(pheno))) model.matrix(~ batch)
                else model.matrix(~ batch + pheno)
    } else {
      design <- if (is.null(unlist(pheno))) model.matrix(~ batch + cluster)
                else model.matrix(~ batch + cluster + pheno)
    }
    batch_loc <- grep("^batch", colnames(design))

    adj_pcs <- as.matrix(adj_ortho(
      if (iter == 1) dat_pcs else adj_pcs,
      design, batch_loc, reduce = TRUE
    ))

    opt_km_adj <- Optimal_Clusters_KMeans(
      adj_pcs[cell_loc, , drop = FALSE],
      max_clusters = ncl_min:ncl_max,
      plot_clusters = FALSE, verbose = FALSE, criterion = "silhouette"
    )
    ncl_adj     <- max(which(opt_km_adj == max(opt_km_adj))) + ncl_min - 1
    sil_score   <- max(opt_km_adj)
    km          <- KMeans_arma(adj_pcs, ncl_adj, n_iter = 100,
                               "random_subset", verbose = FALSE)
    cluster_adj <- as.factor(predict_KMeans(adj_pcs, km))
    bic_km_adj  <- Optimal_Clusters_KMeans(
      adj_pcs, max_clusters = ncl_min:ncl_max,
      plot_clusters = FALSE, verbose = FALSE, criterion = "BIC"
    )

    accepted <- bic_km_adj[ncl_adj - ncl_min + 1] <
                bic_km[ncl_adj - ncl_min + 1]

    convergence_history <- rbind(convergence_history, data.frame(
      iter = iter, k = ncl_adj,
      bic = bic_km_adj[ncl_adj - ncl_min + 1],
      silhouette = sil_score, accepted = accepted
    ))

    if (!accepted) {
      break
    } else {
      cat(paste0("  Iteration ", iter, ": ", ncl_adj, " clusters\n"))
      iter    <- iter + 1
      bic_km  <- bic_km_adj
      cluster <- cluster_adj
    }
  }

  list(cluster = cluster, n_pcs = n_pcs,
       history = convergence_history, elbow_res = elbow_res)
}
