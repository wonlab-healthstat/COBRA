# Internal helper: solve normal equations column-by-column
lm.fit.large <- function(xpp, dat) {
  xppt    <- t(xpp)
  xpptxpp <- as.matrix(xppt %*% xpp)

  A <- if (is.infinite(determinant(xpptxpp)$modulus)) {
    MASS::ginv(xpptxpp) %*% xppt
  } else {
    tryCatch(chol2inv(chol(xpptxpp)) %*% xppt,
             error = function(e) MASS::ginv(xpptxpp) %*% xppt)
  }

  coef <- matrix(NA_real_, nrow = ncol(xpp), ncol = ncol(dat))
  for (i in seq_len(ncol(dat))) {
    coef[, i] <- as.matrix(A %*% dat[, i])
  }
  coef
}

# Internal helper: Rfast variant for large matrices
lm.fit.large.new <- function(xpp, dat) {
  xppt    <- Rfast::transpose(xpp)
  xpptxpp <- xppt %*% xpp

  A <- if (is.infinite(determinant(xpptxpp)$modulus)) {
    MASS::ginv(xpptxpp) %*% xppt
  } else {
    tryCatch(chol2inv(chol(xpptxpp)) %*% xppt,
             error = function(e) MASS::ginv(xpptxpp) %*% xppt)
  }

  coef <- matrix(NA_real_, nrow = ncol(xpp), ncol = ncol(dat))
  for (i in seq_len(ncol(dat))) {
    coef[, i] <- A %*% dat[, i]
  }
  coef
}

# Internal: orthogonal projection for batch effect removal
# Reference: Melo et al. (2009) Comp Stat & Data Analysis
adj_ortho <- function(dat, design, loc, reduce) {
  n   <- nrow(design)
  xp  <- design[, loc]
  xnp <- as.matrix(design[, -loc])

  uniq_loc <- which(apply(xnp, 2, function(x) length(unique(x))) == 1)
  uniq_loc <- uniq_loc[!grepl("*Intercept*", names(uniq_loc))]
  if (length(uniq_loc) > 0) xnp <- xnp[, -uniq_loc]

  xnpt     <- t(xnp)
  xnptxnp  <- xnpt %*% xnp
  xnptxnpi <- chol2inv(chol(xnptxnp))
  xpp      <- as.matrix(
    Diagonal(n) %*% xp - xnp %*% (xnptxnpi %*% xnpt %*% xp)
  )

  coef <- if (reduce) {
    lm.fit(xpp, dat)$coefficients
  } else {
    lm.fit.large.new(xpp, dat)
  }

  na_loc <- which(is.na(rowSums(coef)))

  if (length(na_loc) == 0) {
    dat - xpp %*% coef
  } else {
    dat - xpp[, -na_loc, drop = FALSE] %*% coef[-na_loc, , drop = FALSE]
  }
}
