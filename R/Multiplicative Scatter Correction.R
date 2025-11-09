#' MSC correction within samples (in-place, global assignment)
#'
#' Performs Multiplicative Scatter Correction (MSC) on a full dataset,
#' correcting each scan using the mean spectrum of its own sample as reference.
#' The input matrix `data` is overwritten globally with its corrected values.
#'
#' @param data matrix [n_scans x n_wavelengths]; each row = one scan
#' @param class_vec vector of length n_scans giving the class label for each scan
#' @param scans_per_sample integer; number of scans per sample (e.g. 5)
#'
#' @details
#' For each sample, its mean spectrum is computed and used as the reference.
#' Each scan in that sample is linearly adjusted: (scan - a) / b,
#' where a and b come from lm(scan ~ ref). This removes additive and
#' multiplicative differences between scans within the same sample.
#'
#' @example
#' MSC_within_sample(data, class_vec, 5)
#' # After running, data, sample_means, and class_per_sample exist in global env.
#'
#' @export
MSC <- function(data, class_vec, scans_per_sample) {
  n <- nrow(data)
  m <- ncol(data)
  n_samples <- n / scans_per_sample
  if (n %% scans_per_sample != 0)
    stop("Total scans not divisible by scans_per_sample.")

  sample_idx <- rep(seq_len(n_samples), each = scans_per_sample)
  sample_means <- matrix(NA_real_, n_samples, m)
  class_per_sample <- tapply(class_vec, sample_idx, function(x) x[1])

  for (s in seq_len(n_samples)) {
    rows <- which(sample_idx == s)
    ref <- colMeans(data[rows, , drop = FALSE], na.rm = TRUE)
    sample_means[s, ] <- ref
    for (r in rows) {
      fit <- lm(data[r, ] ~ ref)
      a <- coef(fit)[1]; b <- coef(fit)[2]
      data[r, ] <- (data[r, ] - a) / b
    }
  }

  assign("data", data, envir = .GlobalEnv)
  assign("sample_means", sample_means, envir = .GlobalEnv)
  assign("class_per_sample", class_per_sample, envir = .GlobalEnv)
  invisible()
}
