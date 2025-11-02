#' @title Compute Local Rank
#' @description Estimates the local numerical rank of spectral windows using singular value decomposition (SVD).
#' @param X Numeric matrix or data frame of spectra (rows = samples, columns = wavelengths).
#' @param window_mode Character; "static" or "dynamic" window definition method.
#' @param base_window Integer; base window size for non-peak regions.
#' @param peak_window Integer; window size for peak regions (used only in dynamic mode).
#' @param tol Numeric; singular value ratio threshold for rank determination.
#' @return Integer vector of local ranks for each wavelength position.
#' @export
Compute_Local_Rank <- function(X, window_mode = "static", 
                               base_window = 5, peak_window = 20, 
                               tol = 1e-4) {
  nvar <- ncol(X)
  ranks <- integer(nvar)
  
  # --- Determine window size vector ---
  if (window_mode == "dynamic") {
    if (!exists("peaks", envir = .GlobalEnv))
      stop("Dynamic mode requires 'peaks' vector in global environment")
    peaks <- get("peaks", envir = .GlobalEnv)
    window_n <- rep(base_window, nvar)
    window_n[peaks] <- peak_window
  } else if (window_mode == "static") {
    window_n <- rep(base_window, nvar)
  } else {
    stop("window_mode must be 'static' or 'dynamic'")
  }
  
  # --- Compute local rank for each wavelength ---
  for (i in seq_len(nvar)) {
    half <- floor(window_n[i] / 2)
    low  <- max(1, i - half)
    high <- min(nvar, i + half)
    Xwin <- X[, low:high, drop = FALSE]
    
    s <- svd(Xwin, nu = 0, nv = 0)$d
    ratio <- s / s[1]
    ranks[i] <- max(1, sum(ratio > tol))
  }
  
  return(ranks)
}
