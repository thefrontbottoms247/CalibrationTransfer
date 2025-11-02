#' @title Bin spectral data
#' @description Averages every fixed number of wavelength points into bins.
#' @param data Numeric matrix or data frame of spectra (rows = samples, columns = wavelengths)
#' @param wl Numeric vector of wavelengths corresponding to columns in data
#' @param points_to_bin Integer; number of consecutive points to average per bin
#' @return A list with binned wavelength vector and binned data matrix
#' @export
Binning <- function(data, wl, points_to_bin = 3) {
  if (is.list(data)) data <- as.data.frame(data)
  if (is.data.frame(data)) data <- as.matrix(data)
  if (is.null(dim(data))) data <- matrix(data, nrow = 1)
  storage.mode(data) <- "numeric"
  if (!is.numeric(wl)) wl <- as.numeric(wl)
  
  n_points <- ncol(data)
  n_bins <- floor(n_points / points_to_bin)
  if (n_bins <= 0) stop("Invalid number of bins — check wl length and bin size.")
  
  wl.bin <- numeric(n_bins)
  data.bin <- matrix(NA_real_, nrow = nrow(data), ncol = n_bins)
  
  for (j in seq_len(n_bins)) {
    idx <- ((j - 1) * points_to_bin + 1):((j - 1) * points_to_bin + points_to_bin)
    wl.bin[j] <- mean(wl[idx])
    data.bin[, j] <- rowMeans(data[, idx, drop = FALSE])
  }
  
  list(wl.bin = wl.bin, data.bin = data.bin)
}
