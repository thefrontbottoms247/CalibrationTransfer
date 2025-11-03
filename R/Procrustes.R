#' @title Procrustes_Test
#' @description
#' Perform a Procrustes rotation and permutation test between two PCA or PLS
#' result objects (e.g., parent vs. child instrument models).
#' Returns the protest summary and optional visualization of the Procrustes fit.
#'
#' @param pca1 Object containing a scores matrix (e.g., parent PCA/PLS).
#' @param pca2 Object containing a scores matrix (e.g., child PCA/PLS).
#' @param ncomp Number of components to compare.
#' @param perm Number of permutations for significance testing (default 999).
#' @param plot Logical; if TRUE, plot the Procrustes residuals and rotation.
#' @return List containing the `protest` object, summary, and RV coefficient.
#' @export
#' @examples
#' res <- Procrustes_Test(parent.pca, child.pca, ncomp = 5, plot = TRUE)

Procrustes_Test <- function(pca1, pca2, ncomp = 5, perm = 999, plot = FALSE) {
  if (!requireNamespace("vegan", quietly = TRUE))
    stop("Package 'vegan' is required. Please install it first.")
  
  # --- Input validation ---
  if (is.null(pca1$x) || is.null(pca2$x))
    stop("Both objects must contain an '$x' matrix (e.g., PCA or PLS scores).")
  if (ncol(pca1$x) < ncomp || ncol(pca2$x) < ncomp)
    stop("Both score matrices must have at least ncomp columns.")
  
  X1 <- as.matrix(pca1$x[, 1:ncomp, drop = FALSE])
  X2 <- as.matrix(pca2$x[, 1:ncomp, drop = FALSE])
  
  # --- Procrustes permutation test ---
  pro <- vegan::protest(X1, X2, permutations = perm)
  pro_sum <- summary(pro)
  
  # --- Extract key metrics ---
  rv <- vegan::procrustes(X1, X2, symmetric = TRUE)$ss
  corr <- pro_sum$t0
  pval <- pro_sum$signif
  
  result <- list(
    protest = pro,
    summary = pro_sum,
    RV = rv,
    correlation = corr,
    p.value = pval
  )
  
  # --- Optional visualization ---
  if (plot) {
    vegan::plot(pro, kind = 1, main = paste0(
      "Procrustes Rotation (", ncomp, " components)\n",
      "r = ", format(corr, digits = 3),
      ", p = ", format(pval, digits = 3)
    ))
  }
  
  return(result)
}
