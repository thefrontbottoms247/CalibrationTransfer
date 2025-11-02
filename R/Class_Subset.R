#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export
Class_Subset <- function(bag, classes) {
  # bag: vector of class labels (can be numeric, factor, or character)
  # classes: vector of class indices (or names) to extract
  
  # Ensure 'bag' is treated as factor so classes are well-defined
  bag_factor <- as.factor(bag)
  class_counts <- as.numeric(table(bag_factor))
  class_levels <- levels(bag_factor)
  
  # Support both numeric indices or class names
  if (is.numeric(classes)) {
    if (any(classes < 1 | classes > length(class_levels)))
      stop("Invalid class index in 'classes'.")
    target_levels <- class_levels[classes]
  } else {
    # assume character vector of class names
    if (!all(classes %in% class_levels))
      stop("Some specified classes do not exist in 'bag'.")
    target_levels <- classes
  }
  
  # Return indices corresponding to selected classes
  which(bag_factor %in% target_levels)
}
