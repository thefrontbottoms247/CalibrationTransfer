#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export
Leverage <- function(X, ncomp = 10, frac = 0.15, scale_data = TRUE) {
  if (scale_data) X <- scale(X)
  
  # PCA
  pca <- prcomp(X, scale. = FALSE, center = FALSE)
  T <- pca$x[, 1:min(ncomp, ncol(pca$x)), drop = FALSE]
  
  # Compute leverage (hat values)
  H <- T %*% solve(t(T) %*% T) %*% t(T)
  lev <- diag(H)
  
  # Rank samples by leverage
  k <- round(frac * nrow(X))
  subset_idx <- order(lev, decreasing = TRUE)[1:k]
  
  cat(sprintf("Selected %d samples (%.1f%%) with highest leverage.\n", k, frac*100))
  return(sort(subset_idx))
}
