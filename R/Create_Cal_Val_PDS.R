#' @title Create Calibration and Validation Sets for PDS
#' @description
#' Generates calibration and validation subsets for PDS, with optional CORAL correction
#' and feature selection.
#'
#' @param parent Reference spectra (matrix or data.frame)
#' @param child  Target spectra (matrix or data.frame)
#' @param class.vector Class labels corresponding to rows
#' @param wl.vector Wavelength vector
#' @param seed Random seed for reproducibility
#' @param feat_select Optional integer or logical vector for feature selection
#' @param split_method 'partition' (stratified) or 'kfold'
#' @param p Proportion for calibration (if split_method='partition')
#' @param k Number of folds (if split_method='kfold')
#' @param remove_classes Vector of classes to exclude
#' @param ncomp_svd Optional vector of SVD components
#' @param CORAL Logical; if TRUE, applies CORAL correction (child → parent)
#' @return Invisible list of calibration and validation subsets
#' @export
Create_Cal_Val_PDS <- function(parent, child, class.vector, wl.vector,
                               seed, feat_select = NA, split_method = 'partition',
                               p = 0.1, k = 5, remove_classes = NULL,
                               ncomp_svd = NULL, CORAL = FALSE) {
  if (!is.matrix(parent) && !is.data.frame(parent))
    stop("parent must be a matrix or data.frame")
  if (!is.matrix(child) && !is.data.frame(child))
    stop("child must be a matrix or data.frame")
  if (nrow(parent) != nrow(child))
    stop("parent and child must have the same number of rows")
  if (length(class.vector) != nrow(parent))
    stop("class.vector length must match number of rows")

  set.seed(seed)

  if (!is.null(remove_classes) && length(remove_classes) > 0) {
    keep_idx <- !class.vector %in% remove_classes
    parent <- parent[keep_idx, , drop = FALSE]
    child <- child[keep_idx, , drop = FALSE]
    class.vector <- class.vector[keep_idx]
  }

  # --- Feature selection before splitting and CORAL ---
  if (!is.null(feat_select) && !anyNA(feat_select) && is.vector(feat_select)) {
    wl <- wl.vector[feat_select]
    parent <- parent[, feat_select, drop = FALSE]
    child  <- child[, feat_select, drop = FALSE]

    if (!is.null(ncomp_svd)) {
      if (length(ncomp_svd) == length(wl.vector)) {
        ncomp_svd <- ncomp_svd[feat_select]
        cat(sprintf("Note: ncomp_svd vector subsetted from %d to %d variables\n",
                    length(wl.vector), length(ncomp_svd)))
      } else {
        stop("ncomp_svd length mismatch with wl.vector length.")
      }
    }

    if (exists("ncomp_vec", envir = .GlobalEnv)) {
      ncomp_vec_original <- get("ncomp_vec", envir = .GlobalEnv)
      if (length(ncomp_vec_original) == length(wl.vector)) {
        ncomp_vec_subset <- ncomp_vec_original[feat_select]
        assign("ncomp_vec", ncomp_vec_subset, envir = .GlobalEnv)
        cat(sprintf("Note: ncomp_vec vector subsetted from %d to %d variables\n",
                    length(wl.vector), length(ncomp_vec_subset)))
      } else {
        warning("ncomp_vec length mismatch – not subsetting.")
      }
    }
  } else {
    wl <- wl.vector
  }

  # --- Optional CORAL correction (applied after feature selection) ---
  if (CORAL) {
    cat("Applying CORAL alignment (child → parent)...\n")
    child <- CORAL(child, parent,
                   center = TRUE,
                   scale = FALSE,
                   lambda = 1e-6,
                   method = "svd",
                   rank = NULL)
  }

  # --- Calibration / Validation split ---
  if (split_method == "partition") {
    if (is.null(p)) stop("For 'partition', specify p (0–1).")
    ind_cal <- createDataPartition(class.vector, p = p, list = FALSE)
    ind_val <- setdiff(seq_along(class.vector), ind_cal)
  } else if (split_method == "kfold") {
    if (is.null(k)) stop("For 'kfold', specify k.")
    folds <- createFolds(class.vector, k = k, list = FALSE)
    ind_cal <- which(folds != 1)
    ind_val <- which(folds == 1)
  } else stop("split_method must be 'partition' or 'kfold'.")

  parent.trans <- as.matrix(parent[ind_cal, , drop = FALSE])
  child.trans  <- as.matrix(child[ind_cal, , drop = FALSE])
  parent.mod   <- as.matrix(parent[ind_val, , drop = FALSE])
  child.mod    <- as.matrix(child[ind_val, , drop = FALSE])

  class_cal <- class.vector[ind_cal]
  class_val <- class.vector[ind_val]

  parent.trans <- data.frame(parent.trans); colnames(parent.trans) <- wl
  child.trans  <- data.frame(child.trans);  colnames(child.trans)  <- wl
  parent.mod   <- data.frame(parent.mod);   colnames(parent.mod)   <- wl
  child.mod    <- data.frame(child.mod);    colnames(child.mod)    <- wl

  parent.mod$Class <- as.factor(class_val)
  child.mod$Class  <- as.factor(class_val)

  assign("parent.trans", parent.trans, envir = .GlobalEnv)
  assign("child.trans",  child.trans,  envir = .GlobalEnv)
  assign("parent.mod",   parent.mod,   envir = .GlobalEnv)
  assign("child.mod",    child.mod,    envir = .GlobalEnv)
  assign("class_cal",    class_cal,    envir = .GlobalEnv)
  assign("class_val",    class_val,    envir = .GlobalEnv)
  assign("ind_PDS_trans", ind_cal, envir = .GlobalEnv)
  assign("ind_PDS_model", ind_val, envir = .GlobalEnv)

  if (!is.null(ncomp_svd))
    assign("ncomp_svd", ncomp_svd, envir = .GlobalEnv)

  invisible(list(
    parent.trans = parent.trans,
    child.trans  = child.trans,
    parent.mod   = parent.mod,
    child.mod    = child.mod,
    wl = wl,
    class_cal = class_cal,
    class_val = class_val,
    ncomp_svd = ncomp_svd
  ))
}
