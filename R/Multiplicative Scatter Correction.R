#' @title Multiplicative Scatter Correction (MSC)
#' @description
#' Performs multiplicative scatter correction (MSC) on spectral data to reduce
#' additive and multiplicative effects caused by light scattering and particle size
#' variations between samples or instruments.
#'
#' @param X A numeric matrix or data frame of sample spectra (rows = samples, columns = variables).
#' @param reference Optional numeric vector or matrix representing the reference spectrum.
#'                  If not provided, the mean spectrum of \code{X} is used.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Xmsc}}{MSC-corrected spectra.}
#'   \item{\code{reference}}{Reference spectrum used for correction.}
#'   \item{\code{coefficients}}{Matrix of regression intercepts and slopes (columns: intercept, slope).}
#' }
#'
#' @details
#' For each sample spectrum \eqn{x_i}, MSC performs a simple linear regression against
#' the reference spectrum \eqn{x_{ref}}:
#' \deqn{x_i = a_i + b_i x_{ref} + e_i}
#' and corrects the spectrum as:
#' \deqn{x_{corr,i} = (x_i - a_i) / b_i}
#'
#' @examples
#' \dontrun{
#' result <- MSC(X = spectra_matrix)
#' Xcorr <- result$Xmsc
#' }
#'
#' @export
MSC <- function(X, reference = NULL) {
  X <- as.matrix(X)
  
  # --- Define reference spectrum ---
  if (is.null(reference)) {
    reference <- colMeans(X, na.rm = TRUE)
  } else {
    reference <- as.numeric(reference)
    if (length(reference) != ncol(X))
      stop("Reference spectrum length must match number of columns in X.")
  }
  
  # --- Initialize result matrices ---
  Xmsc <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
  coefs <- matrix(NA_real_, nrow = nrow(X), ncol = 2,
                  dimnames = list(NULL, c("Intercept", "Slope")))
  
  # --- Apply MSC sample-wise ---
  for (i in seq_len(nrow(X))) {
    fit <- stats::lm(X[i, ] ~ reference)
    a <- coef(fit)[1]
    b <- coef(fit)[2]
    coefs[i, ] <- c(a, b)
    Xmsc[i, ] <- (X[i, ] - a) / b
  }
  
  list(
    Xmsc = Xmsc,
    reference = reference,
    coefficients = coefs
  )
}
