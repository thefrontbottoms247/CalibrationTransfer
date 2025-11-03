#' @title Orthogonal_Procrustes
#' @description
#' Compute the optimal orthogonal (rotation + uniform scaling + translation)
#' transformation that aligns matrix Y to matrix X using the classic
#' Orthogonal Procrustes method (minimizing Frobenius norm).
#'
#' @param X Reference matrix (n × p).
#' @param Y Target matrix to align (n × p).
#' @return
#' A list with:
#' \item{R}{Rotation matrix (p × p).}
#' \item{scale}{Uniform scaling factor.}
#' \item{trans}{Translation vector.}
#' \item{transform}{Closure function for transforming new data.}
#' @export
#' @examples
#' op <- Orthogonal_Procrustes(Xref, Xchild)
#' X_aligned <- op$transform(Xchild)

Orthogonal_Procrustes <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (!all(dim(X) == dim(Y)))
    stop("X and Y must have identical dimensions.")
  if (nrow(X) < ncol(X))
    warning("Fewer observations than dimensions; results may be unstable.")
  
  # --- Center each dataset ---
  cx <- colMeans(X)
  cy <- colMeans(Y)
  Xc <- sweep(X, 2, cx)
  Yc <- sweep(Y, 2, cy)
  
  # --- Core SVD step ---
  S <- t(Xc) %*% Yc
  sv <- svd(S)
  R <- sv$u %*% t(sv$v)
  
  # --- Correct improper rotation (reflection) ---
  if (det(R) < 0) {
    sv$u[, ncol(sv$u)] <- -sv$u[, ncol(sv$u)]
    R <- sv$u %*% t(sv$v)
  }
  
  # --- Scale & translation ---
  s <- sum(sv$d) / sum(Xc^2)
  tvec <- as.numeric(cy - s * (cx %*% R))
  
  # --- Define reusable transform closure ---
  transform <- function(Xnew) {
    Xnew <- as.matrix(Xnew)
    s * (Xnew %*% R) + matrix(tvec, nrow(Xnew), ncol(Xnew), byrow = TRUE)
  }
  
  structure(list(
    R = R,
    scale = s,
    trans = tvec,
    transform = transform
  ), class = "OrthogonalProcrustes")
}
