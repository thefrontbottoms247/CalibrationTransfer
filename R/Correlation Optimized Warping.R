#' @title Correlation Optimized Warping (COW)
#' @description
#' Performs correlation optimized warping (COW) for spectral alignment.
#' The algorithm divides spectra into segments and warps them within a
#' specified slack to maximize correlation with a reference spectrum.
#'
#' @param X A numeric matrix or data frame of spectra (rows = samples, columns = variables).
#' @param reference Optional numeric vector or matrix specifying the reference spectrum.
#'                  If NULL, the mean spectrum of \code{X} is used.
#' @param segment Integer. Number of segments to divide each spectrum into.
#' @param slack Integer. Maximum number of variables by which each segment can be stretched or compressed.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Xcow}}{COW-aligned spectra.}
#'   \item{\code{reference}}{Reference spectrum used for alignment.}
#'   \item{\code{warping}}{Warping paths for each sample.}
#' }
#'
#' @details
#' This function uses the implementation from the \pkg{prospectr} package.
#' COW is typically used before building PLS or transfer models to correct
#' wavelength shifts and improve model transferability.
#'
#' @examples
#' \dontrun{
#' result <- COW(parent_spectra, segment = 20, slack = 2)
#' aligned <- result$Xcow
#' }
#'
#' @export
COW <- function(X, reference = NULL, segment = 20, slack = 2) {
  if (!requireNamespace("prospectr", quietly = TRUE))
    stop("Package 'prospectr' required but not installed.")
  
  X <- as.matrix(X)
  if (is.null(reference))
    reference <- colMeans(X, na.rm = TRUE)
  
  cow_result <- prospectr::cow(Xr = reference, X = X, segment = segment, slack = slack)
  
  list(
    Xcow = cow_result$Xwarped,
    reference = reference,
    warping = cow_result$warp.info
  )
}
