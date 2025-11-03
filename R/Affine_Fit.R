#' @title Affine_Fit
#' @description
#' Compute the optimal affine transformation (linear + translation)
#' that maps matrix X onto matrix Y by least squares minimization:
#' \deqn{ \min_{B,t} \| X B + \mathbf{1}t^\top - Y \|_F }.
#' This allows non-orthogonal scaling, shear, and translation.
#'
#' @param X Reference matrix (n × p).
#' @param Y Target matrix (n × p).
#' @return
#' A list with:
#' \item{B}{Coefficient matrix including translation (p + 1 × p).  
#'          The last row corresponds to translation.}
#' \item{A}{Linear transformation submatrix (p × p).}
#' \item{tvec}{Translation vector (1 × p).}
#' \item{transform}{Closure function for transforming new data.}
#' @export
#' @examples
#' af <- Affine_Fit(Xref, Xchild)
#' X_aligned <- af$transform(Xchild)

Affine_Fit <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # --- Sanity checks ---
  if (!all(dim(X) == dim(Y)))
    stop("X and Y must have identical dimensions.")
  if (nrow(X) <= ncol(X))
    warning("Fewer observations than dimensions; affine solution may be unstable.")
  
  # --- Compute affine coefficients (least squares) ---
  B <- qr.solve(cbind(X, 1), Y)  # (p+1) × p matrix
  
  # --- Split into linear & translation parts ---
  A <- B[1:ncol(X), , drop = FALSE]
  tvec <- B[ncol(X) + 1, , drop = FALSE]
  
  # --- Reusable transformation function ---
  transform <- function(Xnew) {
    Xnew <- as.matrix(Xnew)
    cbind(Xnew, 1) %*% B
  }
  
  structure(list(
    B = B,
    A = A,
    tvec = tvec,
    transform = transform
  ), class = "AffineFit")
}
