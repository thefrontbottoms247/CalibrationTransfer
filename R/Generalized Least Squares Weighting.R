#' @title Generalized Least Squares Weighting (GLSW)
#' @description
#' Performs generalized least squares weighting between spectra from
#' source and target instruments to compensate for inter-instrumental variance.
#' This method computes a weighting matrix based on the covariance structure
#' of the source and target datasets and applies it to correct the source spectra.
#'
#' @param X_S A numeric matrix or data frame of source (reference) spectra.
#' @param Y_S A numeric vector of reference concentrations corresponding to \code{X_S}.
#' @param X_T A numeric matrix or data frame of target (secondary) spectra.
#' @param ncomp Integer. Number of PLS components to use (default = 5).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{W}}{Weighting matrix derived from covariance ratio.}
#'   \item{\code{pls_model}}{Trained PLS model on weighted source data.}
#'   \item{\code{Y_pred}}{Predictions on weighted target data.}
#'   \item{\code{X_S_weighted}}{Weighted source spectra.}
#'   \item{\code{X_T_weighted}}{Weighted target spectra.}
#' }
#' @details
#' The weighting matrix is computed as:
#' \deqn{W = cov(X_S)^{-1} cov(X_T)}
#' where \eqn{cov(X_S)} and \eqn{cov(X_T)} are covariance matrices
#' of source and target spectra, respectively. The generalized inverse
#' (\code{MASS::ginv}) is used to ensure numerical stability.
#'
#' @examples
#' \dontrun{
#' result <- GLSW(X_S = parent_spectra,
#'                Y_S = parent_conc,
#'                X_T = child_spectra,
#'                ncomp = 5)
#' result$Y_pred
#' }
#'
#' @export
GLSW <- function(X_S, Y_S, X_T, ncomp = 5) {
  # Dependency checks
  if (!requireNamespace("pls", quietly = TRUE))
    stop("Package 'pls' required but not installed.")
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' required but not installed.")
  
  # --- Input validation ---
  X_S <- as.matrix(X_S)
  X_T <- as.matrix(X_T)
  Y_S <- as.numeric(Y_S)
  
  if (nrow(X_S) != length(Y_S))
    stop("Number of rows in X_S must match length of Y_S.")
  
  # --- Compute covariance-based weighting matrix ---
  cov_S <- stats::cov(X_S)
  cov_T <- stats::cov(X_T)
  W <- MASS::ginv(cov_S) %*% cov_T
  
  # --- Apply weights ---
  X_S_weighted <- X_S %*% W
  X_T_weighted <- X_T %*% W
  
  # --- Train PLS model ---
  pls_model <- pls::plsr(Y_S ~ X_S_weighted, ncomp = ncomp, validation = "CV")
  
  # --- Predict on target data ---
  Y_pred <- predict(pls_model, newdata = X_T_weighted)
  
  # --- Return structured output ---
  result <- list(
    W = W,
    pls_model = pls_model,
    Y_pred = Y_pred,
    X_S_weighted = X_S_weighted,
    X_T_weighted = X_T_weighted
  )
  
  invisible(result)
}
