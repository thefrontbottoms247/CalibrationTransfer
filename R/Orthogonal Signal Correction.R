#' @title Orthogonal Signal Correction (OSC)
#' @description
#' Removes systematic variation in spectral data that is orthogonal
#' (i.e., unrelated) to the response variable.  
#' The method projects out components of \code{X} that do not contribute
#' to explaining \code{Y}, improving model robustness across instruments.
#'
#' @param X A numeric matrix or data frame of spectra (rows = samples, columns = variables).
#' @param Y A numeric vector or matrix of response values corresponding to \code{X}.
#' @param ncomp Integer. Number of orthogonal components to remove (default = 1).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Xcorr}}{OSC-corrected predictor matrix.}
#'   \item{\code{P}}{Matrix of orthogonal loadings.}
#'   \item{\code{T}}{Matrix of orthogonal scores.}
#' }
#'
#' @details
#' Orthogonal Signal Correction removes latent components of \code{X}
#' that are orthogonal to \code{Y}.  
#' At each iteration, an orthogonal component \eqn{t_o} is computed as:
#' \deqn{t_o = X w_o, \quad w_o = (I - P_y P_y^\top) X^\top X w_o}
#' where \eqn{P_y} is the projection of \code{Y}.  
#' The corrected matrix is then:
#' \deqn{X_{corr} = X - t_o p_o^\top}
#'
#' This implementation follows the standard approach used in chemometrics literature
#' (Wold et al., 1998).
#'
#' @references
#' Wold, S., Antti, H., Lindgren, F., & Öhman, J. (1998). Orthogonal signal correction of near-infrared spectra. *Chemometrics and Intelligent Laboratory Systems*, 44(1–2), 175–185.
#'
#' @examples
#' \dontrun{
#' result <- OSC(X = spectra, Y = response, ncomp = 2)
#' Xcorr <- result$Xcorr
#' }
#'
#' @export
OSC <- function(X, Y, ncomp = 1) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  if (nrow(X) != nrow(Y))
    stop("X and Y must have the same number of rows (samples).")
  
  # Center X and Y
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  
  n <- nrow(X)
  p <- ncol(X)
  
  T_osc <- matrix(0, n, ncomp)
  P_osc <- matrix(0, p, ncomp)
  Xcorr <- X
  
  for (a in seq_len(ncomp)) {
    # Compute weight vector orthogonal to Y
    w <- crossprod(Xcorr, Y)
    w <- w / sqrt(sum(w^2))
    
    # Orthogonal component scores
    t <- Xcorr %*% w
    p <- crossprod(Xcorr, t) / drop(crossprod(t))
    p <- p / sqrt(sum(p^2))
    
    # Make component orthogonal to Y
    t_o <- t - (Y %*% (crossprod(Y, t) / drop(crossprod(Y))))
    p_o <- crossprod(Xcorr, t_o) / drop(crossprod(t_o))
    
    # Deflate X
    Xcorr <- Xcorr - t_o %*% t(p_o)
    
    # Store
    T_osc[, a] <- t_o
    P_osc[, a] <- p_o
  }
  
  list(
    Xcorr = Xcorr,
    T = T_osc,
    P = P_osc
  )
}
