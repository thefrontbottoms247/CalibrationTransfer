#' @title External Parameter Orthogonalization (EPO)
#' @description
#' Removes systematic variations between master and slave instruments
#' by projecting their spectral differences onto a reduced orthogonal subspace.
#' The method identifies latent components of the difference matrix and subtracts
#' their contribution from the slave spectra.
#'
#' @param master A numeric matrix or data frame representing the reference (master) spectra.
#' @param slave  A numeric matrix or data frame representing the target (slave) spectra.
#'               Must have the same number of samples (rows) as \code{master}.
#' @param ncomp  Integer. The number of principal components of the difference matrix
#'               to use for correction.
#'
#' @return A corrected slave matrix with orthogonalized signal components removed.
#' @details
#' This implementation follows the standard EPO approach:
#' \deqn{X_{corr} = X_{slave} - T_p T_p^\top}
#' where \eqn{T_p} are the retained principal component scores of the difference matrix.
#'
#' The global variable \code{slave_cor} is assigned for convenience,
#' but the corrected spectra are also returned invisibly for use in pipelines.
#'
#' @examples
#' \dontrun{
#' corrected <- EPO(master_spectra, slave_spectra, ncomp = 5)
#' }
#'
#' @export
EPO <- function(master, slave, ncomp) {
  # Load required packages (optional — avoid inside packages, see note)
  if (!requireNamespace("pracma", quietly = TRUE))
    stop("Package 'pracma' required but not installed.")
  if (!requireNamespace("prospectr", quietly = TRUE))
    stop("Package 'prospectr' required but not installed.")
  if (!requireNamespace("pls", quietly = TRUE))
    stop("Package 'pls' required but not installed.")
  
  master <- as.matrix(master)
  slave  <- as.matrix(slave)
  
  if (nrow(master) != nrow(slave))
    stop("Master and slave matrices must have the same number of samples (rows).")
  
  # --- Compute difference and its PCA ---
  dif_mtrx <- slave - master
  pca_dif  <- prcomp(dif_mtrx, center = TRUE, scale. = FALSE)
  
  # --- Reconstruct systematic variation ---
  cor_mtrx <- pca_dif$x[, 1:ncomp] %*% t(pca_dif$rotation[, 1:ncomp])
  
  # --- Apply correction ---
  slave_cor <- slave - cor_mtrx
  
  assign("slave_cor", slave_cor, envir = .GlobalEnv)  # for session reuse
  
  invisible(slave_cor)
}
