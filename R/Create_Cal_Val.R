#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export
Create_Cal_Val <- function(data.matrix, class.vector, wavelength.vector, 
                           seed, feature_selection, split_method){

  ### reference wood chip sampling
  set.seed(seed)
  
  if (split_method == "partition") {
    p <- as.numeric(readline("Enter partition percentage: "))
    if (is.na(p) || p <= 0 || p >= 1)
      stop("Invalid value for p. Must be a number between 0 and 1.")
    ind <- as.vector(createDataPartition(class.vector, p = p, list = FALSE))
    #### remove as vector when making this have more than 1 partition
    
  } else if (split_method == "kfold") {
    if (is.null(k)) stop("For 'kfold' split_method, you must specify number of folds 'k'.")
    ind <- createFolds(class.vector, k = k, list = FALSE)
  }
  
  cal <- as.matrix(data.matrix[ind, , drop = FALSE])
  val <- as.matrix(data.matrix[-ind, , drop = FALSE])
  
  if (!missing(feature_selection) && is.vector(feature_selection)) {
    wl <- wavelength.vector[feature_selection]
    cal <- cal[,feature_selection]
    val <- val[,feature_selection]
    
    cal <- data.frame(cal)
    cal$Class <- as.factor(class.vector[ind])
    colnames(cal) <- c(wl, "Class")
    
    val <- data.frame(val)
    val$Class <- as.factor(class.vector[-ind])
    colnames(val) <- c(wl, "Class")
   
    assign("cal", cal, envir = .GlobalEnv)
    assign("val", val, envir = .GlobalEnv)
    
    return()
    
  }
    
  cal <- data.frame(cal)
  cal$Class <- as.factor(class.vector[ind])
  colnames(cal) <- c(wavelength.vector, "Class")
  
  val <- data.frame(val)
  val$Class <- as.factor(class.vector[-ind])
  colnames(val) <- c(wavelength.vector, "Class")

  assign("cal", cal, envir = .GlobalEnv)
  assign("val", val, envir = .GlobalEnv)
  
}
