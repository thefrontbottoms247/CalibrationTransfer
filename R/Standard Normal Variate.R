#' @title Standard Normal Variate (SNV) Correction
#' @description
#' Performs Standard Normal Variate (SNV) normalization on both master and slave
#' spectral datasets to remove additive and multiplicative scatter effects.
#'
#' @param master A numeric matrix or data frame of spectra from the master instrument
#'               (rows = samples, columns = wavelengths or variables).
#' @param slave  A numeric matrix or data frame of spectra from the slave instrument,
#'               with the same number of columns as \code{master}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{master_snv}}{SNV-corrected master spectra.}
#'   \item{\code{slave_snv}}{SNV-corrected slave spectra.}
#'   \item{\code{mean_cor}}{Row-wise mean differences between slave and master after SNV.}
#'   \item{\code{scale_cor}}{Row-wise scale ratios between slave and master after SNV.}
#' }
#'
#' @details
#' For each sample (row), SNV standardizes by subtracting its mean and dividing by
#' its standard deviation:
#' \deqn{x_{snv} = (x - \bar{x}) / s_x}
#' This function applies the same operation to both master and slave data and
#' computes the corresponding mean and scale corrections for later use in calibration transfer.
#'
#' @examples
#' \dontrun{
#' result <- SNV(master_spectra, slave_spectra)
#' corrected_master <- result$master_snv
#' }
#'
#' @export
SNV <- function(master, slave) {
  master <- as.matrix(master)
  slave  <- as.matrix(slave)
  
  if (ncol(master) != ncol(slave))
    stop("Master and slave must have the same number of variables (columns).")
  
  # --- Apply SNV normalization row-wise ---
  master_snv <- t(apply(master, 1, function(x) (x - mean(x)) / sd(x)))
  slave_snv  <- t(apply(slave, 1,  function(x) (x - mean(x)) / sd(x)))
  
  # --- Compute correction terms ---
  mean_cor  <- rowMeans(slave_snv) - rowMeans(master_snv)
  scale_cor <- apply(slave_snv, 1, sd) / apply(master_snv, 1, sd)
  
  list(
    master_snv = master_snv,
    slave_snv  = slave_snv,
    mean_cor   = mean_cor,
    scale_cor  = scale_cor
  )
}
