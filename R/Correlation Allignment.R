#' @title CORAL: Correlation Alignment for Domain Adaptation
#' @description
#' Aligns covariance structures between a parent (source) and child (target)
#' instrument for unsupervised calibration transfer. Two computational modes are
#' supported:
#' - **"full"**: classic implementation using full covariance matrices.
#' - **"svd"**: low-rank SVD version for high-dimensional spectra (n << p).
#'
#' @details
#' Given source data \eqn{X_s} and target data \eqn{X_t}, CORAL computes a
#' linear transform \eqn{A = \Sigma_t^{-1/2} \Sigma_s^{1/2}} so that
#' \eqn{X_t A} has the same covariance as \eqn{X_s}.
#'
#' The SVD version operates in the r-dimensional sample subspace (r ≤ n − 1),
#' avoiding the p×p covariance eigendecomposition and reducing memory load from
#' O(p²) to O(p·r).
#'
#' @param child Numeric matrix or data frame from the child (target) instrument.
#' @param parent Numeric matrix or data frame from the parent (source) instrument.
#' @param center Logical; center both datasets before alignment (default TRUE).
#' @param scale Logical; scale variables to unit variance before computing covariances (default FALSE).
#' @param lambda Numeric ridge term added to diagonals for stability (default 1e-6).
#' @param method Character, either `"full"` or `"svd"` (default `"auto"` chooses `"svd"` when n << p).
#' @param rank Optional integer for `"svd"` method, number of retained components.
#'
#' @return Invisibly returns the aligned child spectra with attribute `"A"`
#' (transformation matrix). The aligned data are also assigned globally as `child`.
#'
#' @examples
#' \dontrun{
#' CORAL(child.trans, parent.trans, method = "svd")
#' CORAL(child.trans, parent.trans, method = "full")
#' }
#' @references
#' Sun & Saenko (2016). Deep CORAL: Correlation Alignment for Deep Domain Adaptation.
#' ECCV Workshops, 443–450. Springer.
#'
#' @importFrom MASS ginv
#' @importFrom expm sqrtm
#' @export
CORAL <- function(child, parent,
                  center = TRUE, scale = FALSE,
                  lambda = 1e-6,
                  method = c("auto", "full", "svd"),
                  rank = NULL) {
  method <- match.arg(method)
  child  <- as.matrix(child)
  parent <- as.matrix(parent)

  n_c <- nrow(child); n_p <- nrow(parent)
  p   <- ncol(child)
  if (ncol(parent) != p) stop("Child and parent must have same number of variables.")

  # Decide mode automatically
  if (method == "auto") {
    method <- if (p > (n_c + n_p)) "svd" else "full"
  }

  # --- Preprocessing ---
  if (center) {
    mu_c <- colMeans(child); mu_p <- colMeans(parent)
    child  <- sweep(child,  2, mu_c)
    parent <- sweep(parent, 2, mu_p)
  }
  if (scale) {
    sd_c <- apply(child, 2, sd); sd_p <- apply(parent, 2, sd)
    sd_c[sd_c == 0] <- 1; sd_p[sd_p == 0] <- 1
    child  <- sweep(child,  2, sd_c, "/")
    parent <- sweep(parent, 2, sd_p, "/")
  }

  if (method == "full") {
    # ===============================================================
    # Classic full-covariance CORAL
    # ===============================================================
    Sigma_c <- cov(child) + diag(lambda, p)
    Sigma_p <- cov(parent) + diag(lambda, p)

    Sigma_c_half_inv <- MASS::ginv(expm::sqrtm(Sigma_c))
    Sigma_p_half     <- expm::sqrtm(Sigma_p)
    A <- Sigma_c_half_inv %*% Sigma_p_half
    child_aligned <- child %*% A

  } else {
    # ===============================================================
    # SVD-based low-rank CORAL
    # ===============================================================
    # Thin SVDs (economy)
    sv_c <- svd(child, nu = 0)
    sv_p <- svd(parent, nu = 0)
    rmax <- min(length(sv_c$d), length(sv_p$d))
    if (is.null(rank)) rank <- min(rmax,  min(n_c - 1, n_p - 1))
    r <- rank

    Vc <- sv_c$v[, 1:r, drop = FALSE]
    Vp <- sv_p$v[, 1:r, drop = FALSE]
    Sc <- sv_c$d[1:r]; Sp <- sv_p$d[1:r]
    invsqrt_Sc <- 1 / sqrt(pmax(Sc^2, lambda))
    sqrt_Sp    <- sqrt(pmax(Sp^2,  lambda))

    # compute transform A implicitly in r-subspace
    A <- Vc %*% (diag(invsqrt_Sc, r) %*% (t(Vc) %*% Vp) %*% diag(sqrt_Sp, r)) %*% t(Vp)
    child_aligned <- child %*% A
  }

  if (center) child_aligned <- sweep(child_aligned, 2, mu_p, "+")

  attr(child_aligned, "A") <- A
  assign("child", child_aligned, envir = .GlobalEnv)
  invisible(child_aligned)
}
