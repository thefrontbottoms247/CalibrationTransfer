#' @title Find Half-Maximum Width
#' @description Computes the full width and asymmetry of a peak at half-maximum intensity.
#' @param x Numeric vector of signal intensities.
#' @param wl Numeric vector of corresponding wavelengths.
#' @param peak_index Integer index of the peak maximum.
#' @param target Numeric target intensity (typically half the peak height).
#' @return List with `width`, `left`, `right`, and `asymmetry`.
#' @export
Find_Halfmax_Width <- function(x, wl, peak_index, target) {
  left_idx <- peak_index
  while (left_idx > 1 && x[left_idx] > target) left_idx <- left_idx - 1
  left_wl <- approx(c(x[left_idx], x[left_idx + 1]),
                    c(wl[left_idx], wl[left_idx + 1]),
                    xout = target)$y
  right_idx <- peak_index
  while (right_idx < length(x) && x[right_idx] > target) right_idx <- right_idx + 1
  right_wl <- approx(c(x[right_idx - 1], x[right_idx]),
                     c(wl[right_idx - 1], wl[right_idx]),
                     xout = target)$y
  width <- abs(right_wl - left_wl)
  asymmetry <- abs((wl[peak_index] - left_wl) - (right_wl - wl[peak_index]))
  list(width = width, left = left_wl, right = right_wl, asymmetry = asymmetry)
}

#' @title Lorentzian Equation
#' @description Evaluates a Lorentzian curve at given wavelengths based on area, width, and center parameters.
#' @param A Numeric; integrated area under the Lorentzian peak.
#' @param w Numeric; full width at half maximum (FWHM).
#' @param xc Numeric; center wavelength.
#' @param x Numeric vector of wavelengths.
#' @return Numeric vector of Lorentzian intensity values.
#' @export
Lorentzian_Equation <- function(A, w, xc, x) {
  A0 <- 2 * A / (pi * w)
  y <- A0 * (0.5 * w)^2 / ((x - xc)^2 + (0.5 * w)^2)
  return(y)
}

#' @title Lorentzian Peak Parameter Extraction
#' @description Detects peaks and estimates Lorentzian parameters including center, width, area, and asymmetry.
#' @param test Numeric vector of spectral intensities.
#' @param wl Numeric vector of corresponding wavelengths.
#' @param window Optional numeric range for peak search.
#' @param minheight Minimum peak height threshold.
#' @return Data frame containing PeakID, center (`xc`), area (`A`), width (`w`), asymmetry, and bounds.
#' @export
Lorentzian <- function(test, wl, window = NULL, minheight) {
  test <- as.numeric(test)
  peaks <- pracma::findpeaks(test, minpeakheight = minheight)
  if (is.null(peaks) || nrow(peaks) == 0)
    stop("No peaks found above threshold. Try lowering 'minheight'.")
  
  n_peaks <- nrow(peaks)
  results <- data.frame(PeakID = integer(), xc = numeric(), A = numeric(),
                        w = numeric(), Asymmetry = numeric(),
                        Left = numeric(), Right = numeric(), stringsAsFactors = FALSE)
  
  message("Found ", n_peaks, " peaks. Extracting parameters...")
  
  for (i in seq_len(n_peaks)) {
    ind1 <- peaks[i, 3]; ind2 <- peaks[i, 4]
    xc <- wl[peaks[i, 2]]; yc <- peaks[i, 1]
    wl.trunc <- wl[ind1:ind2]; test.trunc <- test[ind1:ind2]
    HM <- yc * 0.5
    
    width_info <- find_halfmax_width(test.trunc, wl.trunc, peaks[i, 2] - ind1 + 1, HM)
    f <- splinefun(wl.trunc, test.trunc, method = "natural")
    A <- integrate(f, lower = min(wl.trunc), upper = max(wl.trunc))$value
    
    results <- rbind(results, data.frame(PeakID = i, xc = xc, A = A,
                                         w = width_info$width,
                                         Asymmetry = width_info$asymmetry,
                                         Left = width_info$left,
                                         Right = width_info$right))
  }
  
  message("Lorentzian parameter extraction complete.")
  return(results)
}
