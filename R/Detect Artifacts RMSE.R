#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export

Detect_Artefacts_RMSE <- function(parent, child.cor, iqr_factor = 3, min_frac = 0.25) {
  # Remove Class column if present and convert to matrix
  if (is.data.frame(parent) && "Class" %in% colnames(parent)) {
    parent <- parent[, -which(colnames(parent) == "Class"), drop = FALSE]
  }
  if (is.data.frame(child.cor) && "Class" %in% colnames(child.cor)) {
    child.cor <- child.cor[, -which(colnames(child.cor) == "Class"), drop = FALSE]
  }
  
  parent <- as.matrix(parent)
  child.cor <- as.matrix(child.cor)
  
  resid_mat <- parent - child.cor  
  q1 <- apply(resid_mat, 2, quantile, 0.25, na.rm = TRUE)
  q3 <- apply(resid_mat, 2, quantile, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - iqr_factor * iqr
  upper <- q3 + iqr_factor * iqr
  artefact_mask <- t(t(resid_mat) < lower | t(resid_mat) > upper)
  frac_exceed <- colMeans(artefact_mask, na.rm = TRUE)
  artefact_cols <- which(frac_exceed > min_frac)
  
  # Assign to global environment
  resid_mat <<- resid_mat
  artefact_mask <<- artefact_mask
  frac_exceed <<- frac_exceed
  artefact_cols <<- artefact_cols
  
  return(artefact_cols)
}
