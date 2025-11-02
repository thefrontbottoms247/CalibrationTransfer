#' @title Kennard–Stone Sampling
#' @description Selects a representative subset of samples using the Kennard–Stone algorithm.
#' @param data Numeric matrix or data frame of samples (rows = samples, columns = variables).
#' @param k Integer; number of samples to select.
#' @return Integer vector of selected row indices.
#' @export
Kennard_Stone <- function(data, k) {
  Xs <- scale(data)
  D <- as.matrix(dist(Xs))
  idx <- which(D == max(D), arr.ind = TRUE)[1, ]
  subset <- idx
  remaining <- setdiff(seq_len(nrow(data)), subset)
  
  while (length(subset) < k) {
    dmin <- apply(D[remaining, subset, drop = FALSE], 1, min)
    next_idx <- remaining[which.max(dmin)]
    subset <- c(subset, next_idx)
    remaining <- setdiff(remaining, next_idx)
  }
  
  subset
}
