#' @title Onion Sampling for Calibration Subset Selection
#' @description
#' Implements the Onion Sampling method for selecting representative samples
#' from a multivariate dataset. Samples are iteratively selected from the outer
#' layers of PCA score space to ensure uniform coverage of the data distribution.
#'
#' @param X A numeric matrix or data frame (rows = samples, columns = variables).
#' @param ncomp Integer. Number of principal components to use (default = 5).
#' @param n_select Integer. Number of samples to select.
#' @param scale Logical. Whether to scale variables before PCA (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{selected_index}}{Indices of selected samples.}
#'   \item{\code{scores}}{PCA scores used for selection.}
#'   \item{\code{layers}}{Layer assignment (1 = outermost).}
#' }
#'
#' @details
#' Onion Sampling selects samples in layers (“onion shells”) from outermost to innermost
#' regions of the PCA score space. It ensures representative calibration subsets
#' that capture both edge and central variability, unlike purely distance-based
#' methods like Kennard–Stone.
#'
#' @examples
#' \dontrun{
#' result <- OnionSampling(X = spectra, ncomp = 5, n_select = 100)
#' cal_set <- spectra[result$selected_index, ]
#' }
#'
#' @export
OnionSampling <- function(X, ncomp = 5, n_select = 100, scale = TRUE) {
  X <- as.matrix(X)
  if (n_select > nrow(X))
    stop("n_select cannot exceed number of samples.")
  
  # --- PCA model ---
  pca <- prcomp(X, center = TRUE, scale. = scale)
  scores <- pca$x[, 1:ncomp, drop = FALSE]
  
  # --- Distance to origin in score space ---
  d <- sqrt(rowSums(scores^2))
  
  selected <- integer(0)
  remaining <- seq_len(nrow(X))
  layer <- integer(length(d))
  layer_counter <- 1
  
  # --- Iteratively peel off layers ---
  while (length(selected) < n_select && length(remaining) > 0) {
    # Compute percentile threshold for outer shell (top 20% by default)
    thresh <- quantile(d[remaining], probs = 0.8, na.rm = TRUE)
    shell <- remaining[d[remaining] >= thresh]
    
    # If fewer than needed, take all outer points
    add_now <- shell
    layer[add_now] <- layer_counter
    
    selected <- c(selected, add_now)
    remaining <- setdiff(remaining, add_now)
    layer_counter <- layer_counter + 1
  }
  
  # Trim to exact target number (if overshot)
  if (length(selected) > n_select)
    selected <- selected[seq_len(n_select)]
  
  list(
    selected_index = selected,
    scores = scores,
    layers = layer
  )
}
