#' @title Generalized Machine Learning Wrapper
#' @description Trains and evaluates multiple supervised learning models with support for cross-validation, repeated runs, and parallel processing. Automatically generates default tuning grids for many common algorithms.
#' @param datamatrix Numeric matrix or data frame of predictors (rows = samples, columns = features).
#' @param class_vector Factor or vector of class labels.
#' @param train_control_method Character; resampling strategy ("cv", "repeatedcv", etc.).
#' @param train_method Character; caret model identifier (e.g., "svmLinear", "rf", "pls", "glmnet").
#' @param n_runs Integer; number of repeated model training iterations.
#' @param grid Optional tuning grid for model parameters. If NULL, a default grid is generated.
#' @param folds Integer; number of cross-validation folds.
#' @param repeats Integer; number of repetitions for repeated cross-validation.
#' @param parallel Logical; whether to enable parallel processing.
#' @return Invisibly returns a list of caret model objects and saves performance metrics (`metrics`, `models`, `accuracy.cv`, `kappa.cv`) and seeds to the global environment.
#' @export
Machine_Learning <- function(datamatrix, 
                             class_vector, 
                             train_control_method = repeatedcv,
                             train_method = 'svmLinear',
                             n_runs = 5,
                             grid = NULL,
                             folds = 10,
                             repeats = 5,
                             parallel = TRUE) {
  
  invisible(lapply(c(
    # Core
    "caret", "doParallel", "parallel",
    # SVM / Kernel
    "kernlab", "krls", "lssvm",
    # Linear / Regularized
    "glmnet", "MASS",
    # Latent Variable / Chemometric
    "pls", "spls",
    # Tree-Based / Ensemble
    "rpart", "randomForest", "partykit", "gbm", "xgboost", "extraTrees",
    # Instance-Based / Local
    "kknn", "kernelknn",
    # Neural / Deep
    "nnet", "RSNNS", "brnn",
    # Discriminant
    "mda",
    # Robust / Miscellaneous
    "earth", "quantreg",
    # Utility
    "Metrics"
  ), function(x) suppressPackageStartupMessages(
    require(x, character.only = TRUE)
  )))
  
  seed.vec <- c()
  
  # --- Normalize aliases for consistency ---
  if (train_method == "randomForest") train_method <- "rf"
  if (train_method == "cppls") train_method <- "pls"
  if (train_method == "svm") train_method <- "svmRadial"
  
  # --- Parallel setup (conditional) ---
  if (parallel) {
    cl <- makeCluster(parallel::detectCores() - 1)
    registerDoParallel(cl)
    message(paste("✓ Parallel processing enabled with", getDoParWorkers(), "cores"))
  } else {
    registerDoSEQ()
    message("✓ Sequential processing enabled (verbose output will be more visible)")
  }
  
  # --- Input cleanup ---
  if (is.matrix(datamatrix)) datamatrix <- as.data.frame(datamatrix)
  class_vector <- as.factor(as.vector(unlist(class_vector)))
  
  # --- Train control setup ---
  if (train_control_method == "repeatedcv") {
    if (is.null(folds)) folds <- as.numeric(readline("How many folds? "))
    if (is.null(repeats)) repeats <- as.numeric(readline("How many repeats? "))
  } else if (train_control_method == "cv") {
    if (is.null(folds)) folds <- as.numeric(readline("How many folds? "))
    repeats <- NULL
  } else {
    folds <- NULL
    repeats <- NULL
  }
  
  # ===============================================================
  # 🔹 Auto-generate default grid
  # ===============================================================
  if (is.null(grid)) {
    grid <- switch(
      train_method,
      
      # --- SVM / Kernel ---
      "svmLinear"     = expand.grid(C = c(0.01, 0.1, 1, 10, 100)),
      "svmLinear2"    = expand.grid(cost = c(0.01, 0.1, 1, 10, 100)),
      "svmPoly"       = expand.grid(degree = 1:5, scale = c(0.001, 0.01, 0.1), C = c(0.1, 1, 10)),
      "svmRadial"     = expand.grid(sigma = c(0.001, 0.01, 0.1), C = c(0.1, 1, 10)),
      "svmRadialCost" = expand.grid(C = c(0.1, 1, 10)),
      "gaussprRadial" = expand.grid(sigma = c(0.001, 0.01, 0.1)),
      "krlsRadial"    = expand.grid(sigma = c(0.001, 0.01, 0.1)),
      "lssvmRadial"   = expand.grid(sigma = c(0.001, 0.01, 0.1), tau = seq(0.1, 1, 0.3)),
      
      # --- Linear / Regularized ---
      "glm"    = NULL,
      "glmnet" = expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-4, 1, 0.5)),
      "ridge"  = expand.grid(lambda = 10^seq(-4, 1, 0.5)),
      "enet"   = expand.grid(fraction = seq(0.1, 1, 0.1)),
      
      # --- Latent Variable / Chemometric ---
      "pls"  = expand.grid(ncomp = 1:30),
      "pcr"  = expand.grid(ncomp = 1:30),
      "spls" = expand.grid(K = 1:10, eta = seq(0.1, 0.9, 0.2)),
      
      # --- Tree-Based & Ensemble ---
      "rpart"    = expand.grid(cp = seq(0.001, 0.1, 0.01)),
      "rf"       = expand.grid(mtry = c(2, 5, 10, 15)),
      "cforest"  = expand.grid(mtry = c(2, 5, 10)),
      "gbm"      = expand.grid(interaction.depth = c(1, 3, 5),
                               n.trees = c(50, 100, 200),
                               shrinkage = 0.1,
                               n.minobsinnode = 10),
      "xgbTree"  = expand.grid(nrounds = c(50, 100, 200),
                               max_depth = c(3, 6),
                               eta = c(0.05, 0.1),
                               gamma = 0,
                               colsample_bytree = 1,
                               min_child_weight = 1,
                               subsample = 1),
      "extraTrees" = expand.grid(mtry = c(2, 5, 10), numRandomCuts = c(1, 3, 5)),
      
      # --- Instance-Based / Local ---
      "knn"       = expand.grid(k = seq(3, 25, by = 2)),
      "kernelknn" = expand.grid(k = seq(3, 25, by = 2), method = c("gaussian", "triangular")),
      
      # --- Neural / Deep ---
      "nnet"  = expand.grid(size = c(1, 3, 5, 7), decay = c(0, 0.1, 0.5)),
      "mlp"   = expand.grid(size = c(3, 5, 7), decay = c(0, 0.1, 0.5)),
      "mlpML" = expand.grid(layer1 = c(3, 5), layer2 = c(3, 5)),
      "brnn"  = NULL,
      
      # --- Discriminant (Polynomial Discriminant Analysis) ---
      "pda" = expand.grid(lambda = seq(0, 0.1, 0.02)),
      
      # --- Robust / Miscellaneous ---
      "earth"    = expand.grid(nprune = seq(2, 30, 4), degree = 1:3),
      "quantreg" = expand.grid(tau = seq(0.1, 0.9, 0.2)),
      
      NULL
    )
  }
  
  if (is.null(grid))
    stop(paste("No tuning grid defined for method:", train_method))
  
  seed.vec <- sample(seq(1,1e6,1),n_runs)
  
  # ===============================================================
  # 🚀 Training loop
  # ===============================================================
  models <- vector("list", n_runs)
  metrics <- matrix(NA, nrow = n_runs, ncol = 5)
  colnames(metrics) <- c("Param", "Acc", "Kap", "AccSD", "KapSD")
  
  for (i in seq_len(n_runs)) {
    message(paste("▶ Running iteration", i, "of", n_runs, "using", train_method, "..."))
    
    trc <- caret::trainControl(method = train_control_method,
                               number = folds,
                               repeats = repeats,
                               returnResamp = "all",
                               savePredictions = "final",
                               verboseIter = TRUE,
                               allowParallel = parallel)
    
    set.seed(seed.vec[i])
    
    model <- caret::train(
      x = datamatrix,
      y = class_vector,
      method = train_method,
      trControl = trc,
      tuneGrid = grid,
      verbose = TRUE
    )
    
    models[[i]] <- model
    best_idx <- which(rownames(model$results) == rownames(model$bestTune))
    res <- model$results[best_idx, c("Accuracy", "Kappa", "AccuracySD", "KappaSD")]
    
    best_params <- paste(names(model$bestTune), model$bestTune, sep = "=", collapse = ", ")
    metrics[i, 1] <- best_params
    metrics[i, 2:5] <- as.numeric(res)
  }
  
  print(metrics)
  
  if (parallel) {
    stopCluster(cl)
    registerDoSEQ()
  }
  
  kappa.cv <- sapply(models, function(x) mean(x$results$Kappa))
  accuracy.cv <- sapply(models, function(x) mean(x$results$Accuracy))
  
  assign("seed.vec", seed.vec, envir = .GlobalEnv)
  assign("kappa.cv", kappa.cv, envir = .GlobalEnv)
  assign("accuracy.cv", accuracy.cv, envir = .GlobalEnv)
  assign("models", models, envir = .GlobalEnv)
  assign("metrics", metrics, envir = .GlobalEnv)
  
  cat("Average Kappa of Cross-Validation:", mean(kappa.cv), "\n")
  message("✅ Metrics matrix saved to global environment as 'metrics'.")
}