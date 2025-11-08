#' @title Direct Standardization (DS)
#' @description Stable DS: computes T such that child %*% T ≈ parent, then applies to child.mod.
#' @param parent numeric matrix/data.frame (n × p)
#' @param child  numeric matrix/data.frame (n × p) — same samples as parent
#' @param child.mod data.frame to correct; last column is Class
#' @param lambda ridge term (default 1e-3)
#' @param intercept logical; add intercept column (default TRUE)
#' @return Invisibly returns T; assigns child.mod_cor globally.
#' @export
DS <- function(parent, child, child.mod, lambda = 1e-3, intercept = TRUE) {
  parent <- as.matrix(parent)
  child  <- as.matrix(child)
  if (nrow(parent) != nrow(child)) stop("Parent and child must have same #rows (samples).")

  X <- if (intercept) cbind(1, child) else child
  pX <- ncol(X)

  XtX <- crossprod(X)
  XtY <- crossprod(X, parent)
  M   <- XtX + diag(lambda, pX)

  if (rcond(M) < 1e-12) {
    message("Matrix ill-conditioned — using SVD-based pseudo-inverse.")
    sv <- svd(M)
    M_inv <- sv$u %*% diag(1 / sv$d) %*% t(sv$u)  # no extra +lambda here
    T <- M_inv %*% XtY
  } else {
    T <- solve(M, XtY)
  }

  # Apply to child.mod (expect last column = Class)
  Xmod <- as.matrix(child.mod[, -ncol(child.mod), drop = FALSE])
  if (intercept) Xmod <- cbind(1, Xmod)

  if (ncol(Xmod) != nrow(T)) {
    stop(sprintf("Non-conformable: ncol(child.mod%s)=%d but nrow(T)=%d.",
                 if (intercept) "+intercept" else "",
                 ncol(Xmod), nrow(T)))
  }

  child.mod_cor <- Xmod %*% T
  child.mod_cor <- as.data.frame(child.mod_cor)

  # Optional: align corrected column names to parent wavelengths if available
  colnames(child.mod_cor) <- colnames(parent)

  child.mod_cor$Class <- child.mod$Class
  assign("child.mod_cor", child.mod_cor, envir = .GlobalEnv)

  invisible(T)
}
