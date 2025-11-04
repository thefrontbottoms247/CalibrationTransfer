# ============================================================
#' @title Run Model Confusion Matrices Safely
#' @description
#' Applies a list of caret models to a dataset (matrix or data.frame),
#' automatically attaches class labels if needed, aligns factor levels,
#' handles prediction failures, and returns confusion matrices for each model.
#'
#' @param models List of trained caret models.
#' @param dataset Either:
#'   (1) A data.frame containing predictors and a 'Class' column, or
#'   (2) A matrix/data.frame of predictors (without class).
#' @param class Optional vector of class labels if dataset lacks 'Class'.
#'
#' @return List of caret confusionMatrix objects (or NA for failed models).
#' @export
# ============================================================
Run_Confusion <- function(models, dataset, class = NULL) {
  if (is.matrix(dataset)) dataset <- as.data.frame(dataset)

  # --- Add class column if provided separately ---
  if (!"Class" %in% colnames(dataset)) {
    if (is.null(class)) stop("If dataset has no 'Class' column, supply a class vector.")
    if (length(class) != nrow(dataset))
      stop("Length of class vector must match number of rows in dataset.")
    dataset$Class <- class
  }

  # --- Define reference levels ---
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

    if (all(is.na(pred))) {
      warning(paste("All predictions NA for", model$method))
      return(NA)
    }

    pred_factor <- factor(as.character(pred), levels = levels_ref)
    true_factor <- factor(as.character(true_class), levels = levels_ref)

    caret::confusionMatrix(pred_factor, true_factor)
  }

  # --- Extract predictors & response ---
  data_x <- dataset[, !(colnames(dataset) %in% "Class"), drop = FALSE]
  data_y <- dataset$Class

  # --- Run across all models ---
  cms_list <- lapply(models, get_cm,
                     data = data_x,
                     true_class = data_y,
                     levels_ref = levels_ref)

  valid <- sum(sapply(cms_list, function(x) !is.na(x)[1]))
  message("✅ Completed confusion matrices: ", valid, " / ", length(models))

  return(cms_list)
}
