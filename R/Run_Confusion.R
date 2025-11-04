# ============================================================
#' @title Run Model Confusion Matrices Safely
#' @description
#' Applies a list of caret models to a given dataset, automatically
#' aligns factor levels, handles prediction failures, and returns
#' confusion matrix objects for each model.
#'
#' @param models List of trained caret models
#' @param dataset Data frame containing predictors and Class column
#'
#' @return List of caret confusionMatrix objects (or NA if failed)
#' @export
# ============================================================
Run_Confusion <- function(models, dataset) {
  if (!"Class" %in% colnames(dataset))
    stop("Dataset must contain a 'Class' column.")
  
  # --- Build levels from all models (universal reference) ---
  levels_ref <- sort(unique(as.character(dataset$Class)))
  
  # --- Safe CM generator ---
  get_cm <- function(model, data, true_class, levels_ref) {
    pred <- tryCatch(
      predict(model, newdata = data),
      error = function(e) {
        warning(paste("Prediction failed for model:", model$method, "-", e$message))
        rep(NA, nrow(data))
      }
    )
    
    # If all predictions are NA → skip
    if (all(is.na(pred))) {
      warning(paste("All predictions NA for", model$method))
      return(NA)
    }
    
    pred_factor <- factor(as.character(pred), levels = levels_ref)
    true_factor <- factor(as.character(true_class), levels = levels_ref)
    
    # Return full caret confusion matrix
    caret::confusionMatrix(pred_factor, true_factor)
  }
  
  # --- Extract predictors & response ---
  data_x <- dataset[, !(colnames(dataset) %in% "Class"), drop = FALSE]
  data_y <- dataset$Class
  
  # --- Apply to all models ---
  cms_list <- lapply(models, get_cm,
                     data = data_x,
                     true_class = data_y,
                     levels_ref = levels_ref)
  
  # --- Summarize result ---
  valid <- sum(sapply(cms_list, function(x) !is.na(x)[1]))
  message("✅ Completed confusion matrices: ", valid, " / ", length(models))
  
  return(cms_list)
}
