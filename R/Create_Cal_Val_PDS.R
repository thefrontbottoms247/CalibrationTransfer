#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export
Create_Cal_Val_PDS <- function(parent, child, class.vector, wl.vector,
                               seed, feat_select = NA, split_method = 'partition',
                               p = 0.1, k = 5, remove_classes = NULL, ncomp_svd = NULL) {
  if (!is.matrix(parent) && !is.data.frame(parent))
    stop("parent must be a matrix or data.frame")
  if (!is.matrix(child) && !is.data.frame(child))
    stop("child must be a matrix or data.frame")
  if (nrow(parent) != nrow(child))
    stop("parent and child must have the same number of rows")
  if (length(class.vector) != nrow(parent))
    stop("class.vector length must match number of rows")
  
  set.seed(seed)
  
  if (!is.null(remove_classes) && length(remove_classes) > 0)
    class.vector <- class.vector[!class.vector %in% remove_classes]
  
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
  
  if (!is.null(feat_select) && !anyNA(feat_select) && is.vector(feat_select)) {
    wl <- wl.vector[feat_select]
    parent.trans <- parent.trans[, feat_select, drop = FALSE]
    child.trans  <- child.trans[, feat_select, drop = FALSE]
    parent.mod   <- parent.mod[, feat_select, drop = FALSE]
    child.mod    <- child.mod[, feat_select, drop = FALSE]
    
    # Subset ncomp_svd vector if provided
    if (!is.null(ncomp_svd)) {
      if (length(ncomp_svd) == length(wl.vector)) {
        ncomp_svd <- ncomp_svd[feat_select]
        cat(sprintf("Note: ncomp_svd vector subsetted from %d to %d variables\n",
                    length(wl.vector), length(ncomp_svd)))
      } else {
        stop("ncomp_svd length (", length(ncomp_svd), 
             ") doesn't match wl.vector length (", length(wl.vector), ")")
      }
    }
    
    # Subset ncomp_vec if it exists in global environment
    if (exists("ncomp_vec", envir = .GlobalEnv)) {
      ncomp_vec_original <- get("ncomp_vec", envir = .GlobalEnv)
      if (length(ncomp_vec_original) == length(wl.vector)) {
        ncomp_vec_subset <- ncomp_vec_original[feat_select]
        assign("ncomp_vec", ncomp_vec_subset, envir = .GlobalEnv)
        cat(sprintf("Note: ncomp_vec vector subsetted from %d to %d variables\n",
                    length(wl.vector), length(ncomp_vec_subset)))
      } else {
        warning("ncomp_vec exists but length (", length(ncomp_vec_original), 
                ") doesn't match wl.vector length (", length(wl.vector), ") - not subsetting")
      }
    }
  } else {
    wl <- wl.vector
  }
  
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
  
  # Assign ncomp_svd to global environment if provided
  if (!is.null(ncomp_svd)) {
    assign("ncomp_svd", ncomp_svd, envir = .GlobalEnv)
  }
  
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