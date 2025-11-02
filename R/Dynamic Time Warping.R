#' @title Dynamic Time Warping (DTW) Alignment
#' @description
#' Performs Dynamic Time Warping (DTW) alignment between spectra and a reference
#' to correct non-linear wavelength shifts. Particularly useful for Raman, NIR,
#' and LIBS data with local distortions.
#'
#' @param X A numeric matrix or data frame of spectra (rows = samples, columns = variables).
#' @param reference Optional numeric vector representing the reference spectrum.
#'                  If NULL, the mean spectrum is used.
#' @param window.type Character. Type of warping constraint window
#'                    (e.g., \code{"none"}, \code{"sakoechiba"}, \code{"slantedband"}).
#' @param window.size Integer. Window size for constrained warping (default = 20).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Xdtw}}{DTW-aligned spectra.}
#'   \item{\code{reference}}{Reference spectrum used.}
#'   \item{\code{warp_paths}}{List of warping paths for each sample.}
#' }
#'
#' @details
#' This function uses the \pkg{dtw} package to compute minimal-distance warping paths
#' between each spectrum and the reference. Aligned spectra are resampled accordingly.
#'
#' @examples
#' \dontrun{
#' result <- DTW(parent_spectra, window.type = "sakoechiba", window.size = 15)
#' aligned <- result$Xdtw
#' }
#'
#' @export
DTW <- function(X, reference = NULL,
                window.type = "sakoechiba", window.size = 20) {
  if (!requireNamespace("dtw", quietly = TRUE))
    stop("Package 'dtw' required but not installed.")
  
  X <- as.matrix(X)
  if (is.null(reference))
    reference <- colMeans(X, na.rm = TRUE)
  
  Xdtw <- matrix(NA_real_, nrow = nrow(X), ncol = length(reference))
  paths <- vector("list", nrow(X))
  
  for (i in seq_len(nrow(X))) {
    alignment <- dtw::dtw(X[i, ], reference,
                          keep = TRUE,
                          step.pattern = dtw::symmetric2,
                          window.type = window.type,
                          window.size = window.size)
    Xdtw[i, ] <- approx(seq_along(alignment$index1),
                        X[i, alignment$index1],
                        xout = seq_len(length(reference)))$y
    paths[[i]] <- alignment$index1
  }
  
  list(
    Xdtw = Xdtw,
    reference = reference,
    warp_paths = paths
  )
}
