#' @title CORAL: Correlation Alignment for Domain Adaptation
#' @description
#' Performs unsupervised domain adaptation between a \emph{child} (target) and
#' \emph{parent} (source) instrument by aligning their covariance structures.
#' The method "re-colors" the target data so that its covariance matches the
#' source domain. This helps stabilize calibration transfer when cross-instrument
#' differences are largely second-order (covariance) distortions.
#'
#' @details
#' Given source data \eqn{X_s} and target data \eqn{X_t}, CORAL finds a linear
#' transform \eqn{A = \Sigma_t^{-1/2} \Sigma_s^{1/2}} such that the transformed
#' target data \eqn{X_t A} has the same covariance as \eqn{X_s}.
#'
#' Optionally, both datasets are centered before alignment and the target is
#' shifted to match the source mean afterward. A small ridge term is added to
#' each covariance matrix for numerical stability.
#'
#' @param child Numeric matrix or data frame from the \emph{child} (target) instrument.
#' @param parent Numeric matrix or data frame from the \emph{parent} (source) instrument.
#' @param center Logical; if TRUE (default), center both datasets before alignment.
#' @param scale Logical; if TRUE, scale variables to unit variance before computing covariances.
#' @param lambda Numeric; ridge regularization term added to diagonal (default = 1e-6).
#'
#' @return Invisibly returns the aligned child spectra (matrix) with attribute \code{"A"}
#' containing the applied transformation matrix. The aligned spectra are also assigned
#' globally as \code{child}.
#'
#' @examples
#' \dontrun{
#' CORAL(child.trans, parent.trans)
#' # child is now globally replaced with its aligned version
#' PDS(parent.trans, child, child.mod, window = "dynamic", model_type = "PLSR")
#' }
#'
#' @references
#' Sun, B., & Saenko, K. (2016). Deep CORAL: Correlation Alignment for Deep Domain Adaptation.
#' *ECCV Workshops*, 443–450. Springer.
#'
#' @importFrom MASS ginv
#' @importFrom expm sqrtm
#' @export
CORAL <- function(child, parent, center = TRUE, scale = FALSE, lambda = 1e-6) {
  # --- Package requirements ---
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' required but not installed.")
  if (!requireNamespace("expm", quietly = TRUE))
    stop("Package 'expm' required but not installed.")
  
  # --- Convert to matrices ---
  child  <- as.matrix(child)
  parent <- as.matrix(parent)
  
  # --- Optionally center and/or scale ---
  if (center) {
    mu_c <- colMeans(child)
    mu_p <- colMeans(parent)
    child  <- sweep(child,  2, mu_c)
    parent <- sweep(parent, 2, mu_p)
  }
  if (scale) {
    sd_c <- apply(child, 2, sd)
    sd_p <- apply(parent, 2, sd)
    sd_c[sd_c == 0] <- 1
    sd_p[sd_p == 0] <- 1
    child  <- sweep(child,  2, sd_c, "/")
    parent <- sweep(parent, 2, sd_p, "/")
  }
  
  # --- Compute covariances ---
  Sigma_c <- cov(child) + diag(lambda, ncol(child))
  Sigma_p <- cov(parent) + diag(lambda, ncol(parent))
  
  # --- Compute whitening/recoloring transform ---
  Sigma_c_half_inv <- MASS::ginv(expm::sqrtm(Sigma_c))
  Sigma_p_half     <- expm::sqrtm(Sigma_p)
  A <- Sigma_c_half_inv %*% Sigma_p_half
  
  # --- Apply alignment ---
  child_aligned <- child %*% A
  
  # --- Shift back to parent mean if centered ---
  if (center) {
    child_aligned <- sweep(child_aligned, 2, mu_p, "+")
  }
  
  # --- Assign globally and return ---
  attr(child_aligned, "A") <- A
  assign("child", child_aligned, envir = .GlobalEnv)
  
  invisible(child_aligned)
}
