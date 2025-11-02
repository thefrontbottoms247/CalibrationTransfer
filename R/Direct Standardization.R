#' @title Direct Standardization (DS)
#' @description
#' Computes the direct standardization transfer matrix between a parent and child
#' instrument or dataset, such that \eqn{child * T ≈ parent}.
#'
#' @param parent A numeric matrix or data frame containing spectra (rows = samples, columns = variables)
#'               from the reference (parent) instrument.
#' @param child  A numeric matrix or data frame containing spectra from the target (child) instrument,
#'               with the same number of samples as \code{parent}.
#'
#' @return The transfer matrix \eqn{T} as a numeric matrix (invisibly).
#' @details
#' This function solves for the transfer matrix \eqn{T} using a least-squares approach
#' via QR decomposition (\code{qr.solve(child, parent)}). The result can be applied to
#' transform new child spectra to the parent space using \code{child_new %*% T}.
#'
#' @examples
#' \dontrun{
#' T <- DS(parent_spectra, child_spectra)
#' child_corrected <- child_spectra %*% T
#' }
#'
#' @export
DS <- function(parent, child) {
  # Coerce to numeric matrices
  parent <- as.matrix(parent)
  child  <- as.matrix(child)
  
  # Verify same samples
  if (nrow(parent) != nrow(child))
    stop("Parent and child matrices must have the same number of samples (rows).")
  
  # Compute transfer matrix: child * T ≈ parent
  T <- qr.solve(child, parent)
  
  # Assign T globally for reuse (optional, mostly for convenience in sessions)
  assign("T", T, envir = .GlobalEnv)
  
  invisible(T)
}
