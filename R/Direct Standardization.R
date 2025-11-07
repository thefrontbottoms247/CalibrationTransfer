#' @title Direct Standardization (DS)
#' @description
#' Computes a stable Direct Standardization (DS) transfer matrix such that
#' \eqn{child %*% T ≈ parent}, then applies the transformation to a
#' corresponding model dataset (\code{child.mod}).
#'
#' @param parent A numeric matrix or data frame (samples × variables) from the parent instrument.
#' @param child  A numeric matrix or data frame (samples × variables) from the child instrument.
#' @param child.mod A data frame of child spectra to be corrected, including a \code{Class} column.
#' @param lambda Regularization term (default = 1e-3) for numerical stability.
#' @param intercept Logical; if TRUE (default), adds intercept column.
#'
#' @return Invisibly returns the transfer matrix \eqn{T}.
#'         The corrected dataset is assigned globally as \code{child.mod_cor}.
#'
#' @export
DS <- function(parent, child, child.mod, lambda = 1e-3, intercept = TRUE) {
  parent <- as.matrix(parent)
  child  <- as.matrix(child)
  if (nrow(parent) != nrow(child))
    stop("Parent and child must have the same number of samples (rows).")

  if (intercept)
    child <- cbind(1, child)

  p <- ncol(child)
  XtX <- crossprod(child)
  XtY <- crossprod(child, parent)
  M <- XtX + diag(lambda, p)

  if (rcond(M) < 1e-12) {
    message("Matrix nearly singular — using SVD-based pseudo-inverse.")
    svd_M <- svd(M)
    M_inv <- svd_M$u %*% diag(1 / (svd_M$d + lambda)) %*% t(svd_M$u)
    T <- M_inv %*% XtY
  } else {
    T <- solve(M, XtY)
  }

  # --- Apply transformation to model dataset ---
  child.mod_cor <- as.matrix(child.mod[,-ncol(child.mod)]) %*% T
  child.mod_cor <- as.data.frame(child.mod_cor)
  child.mod_cor$Class <- child.mod$Class

  # Assign globally for downstream modeling
  assign("child.mod_cor", child.mod_cor, envir = .GlobalEnv)

  invisible(T)
}
