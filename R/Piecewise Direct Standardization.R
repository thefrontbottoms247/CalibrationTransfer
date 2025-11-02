#' @title Piecewise Direct Standardization (PDS)
#' @description Performs calibration transfer between parent and child instruments using regression models.
#' @param parent Numeric matrix or data frame of reference spectra.
#' @param child Numeric matrix or data frame of target spectra.
#' @param window Character; "static" or "dynamic" window definition method.
#' @param ncomp Integer; number of components for latent variable models.
#' @param model_type Character; regression type ("OLS", "Ridge", "Lasso", "PLSR", "PCR", or "SVM").
#' @param degree Integer; polynomial degree for nonlinear feature expansion.
#' @param svm_kernel Character; SVM kernel type ("linear", "polynomial", or "radial").
#' @param base_window Integer; base window size for static or non-peak regions.
#' @param peak_window Integer; window size for dynamic peak regions.
#' @param ncomp_svd Optional numeric vector of component counts per variable.
#' @param nonneg Logical; whether to enforce non-negativity on regression coefficients.
#' @return A list with the sparse coefficient matrix (`P`) and intercept vector (`intercept`).
#' @export
PDS <- function(parent, child, window, ncomp = 20, model_type = 'svmLinear', 
                degree = 1, svm_kernel = "radial", base_window = 5, 
                peak_window = 20, ncomp_svd = NULL, nonneg = FALSE) {
  
  library(progress)
  library(pls)
  library(glmnet)
  library(Matrix)
  library(e1071)  # for SVM
  library(nnls)   # for non-negative least squares
  
  # Validate input matrices
  stopifnot(ncol(parent) == ncol(child))
  nvar <- ncol(parent)
  
  # Determine window sizes for each variable
  if (window == "dynamic") {
    # Check if 'peaks' exists in the global environment
    if (!exists("peaks", envir = .GlobalEnv)) {
      stop("For dynamic window, 'peaks' vector must exist in global environment")
    }
    peaks <- get("peaks", envir = .GlobalEnv)
    
    # Ask user for both window sizes
    if (is.null(base_window)) base_window <- as.numeric(readline("Base window size (non-peak regions)? "))
    if (is.null(peak_window)) peak_window <- as.numeric(readline("Peak window size (peak regions)? "))
    
    # Initialize all windows to base size
    window_n <- rep(base_window, nvar)
    
    # Set peak window size for indices in peaks
    window_n[peaks] <- peak_window
    
    cat(sprintf("Dynamic window: %d variables with window=%d, %d with window=%d\n",
                nvar - length(peaks), base_window,
                length(peaks), peak_window))
    
  } else if (window == "static") {
    if (is.null(base_window)) window_size <- as.numeric(readline("Window size? "))
    else window_size <- base_window
    window_n <- rep(window_size, nvar)
    
  } else {
    stop("Choose either 'static' or 'dynamic' window")
  }
  
  # Initialize coefficient and intercept matrices
  P <- matrix(0, nrow = nvar, ncol = nvar)
  intercept <- numeric(nvar)
  
  # Progress bar
  pb <- progress_bar$new(
    format = "  \033[1;32mPDS building [:bar] :percent :elapsed\033[0m",
    total = nvar, clear = FALSE, width = 70
  )
  
  # Iterate across spectral variables
  for (i in seq_len(nvar)) {
    # Use the window size specific to this variable position
    half <- floor(window_n[i] / 2)
    low  <- max(1, i - half)
    high <- min(nvar, i + half)
    
    Xmat <- as.matrix(child[, low:high, drop = FALSE])
    yvec <- as.numeric(parent[, i])
    
    # Store original dimensions
    n_orig <- ncol(Xmat)
    
    # Determine number of components - use SVD vector if provided
    ncomp_actual <- ncomp
    if (!is.null(ncomp_svd) && model_type %in% c("PLSR", "PCR")) {
      # Validate ncomp_svd vector length
      if (length(ncomp_svd) != nvar) {
        stop("ncomp_svd vector must have length equal to number of variables (", nvar, ")")
      }
      # Use the pre-computed ncomp for this variable position
      ncomp_actual <- min(ncomp_svd[i], ncol(Xmat))
    }
    
    # Add polynomial features if degree > 1 (not for SVM)
    if (degree > 1 && model_type != "SVM") {
      Xmat_poly <- Xmat
      for (deg in 2:degree) {
        Xmat_poly <- cbind(Xmat_poly, Xmat^deg)
      }
      Xmat <- Xmat_poly
    }
    
    # --- Model fitting ---
    if (nonneg && model_type %in% c("OLS", "Ridge", "Lasso")) {
      # Non-negative constraint for OLS, Ridge, and Lasso
      if (model_type == "OLS") {
        # Non-negative least squares
        nnls_fit <- nnls(Xmat, yvec)
        b_nonneg <- coef(nnls_fit)
        
        # Calculate intercept
        intercept[i] <- mean(yvec - Xmat %*% b_nonneg)
        P[low:high, i] <- b_nonneg[1:n_orig]
        
      } else {
        # Non-negative Ridge or Lasso using glmnet
        alpha <- if (model_type == "Ridge") 0 else 1
        
        cvfit <- cv.glmnet(Xmat, yvec,
                           alpha = alpha,
                           intercept = TRUE,
                           standardize = FALSE,
                           lower.limits = 0)  # Non-negativity constraint
        b <- as.numeric(coef(cvfit, s = "lambda.min"))
        intercept[i] <- b[1]
        P[low:high, i] <- b[2:(n_orig + 1)]
      }
      
    } else if (model_type %in% c("PLSR", "PCR")) {
      df <- data.frame(y = yvec, Xmat)
      colnames(df) <- c("y", paste0("X", seq_len(ncol(Xmat))))
      
      # Use ncomp_actual (SVD-derived or original)
      ncomp_use <- min(ncomp_actual, ncol(Xmat))
      
      fit <- if (model_type == "PLSR") {
        plsr(y ~ ., data = df,
             ncomp = ncomp_use,
             validation = "none")
      } else {
        pcr(y ~ ., data = df,
            ncomp = ncomp_use,
            validation = "none")
      }
      
      b <- drop(coef(fit, ncomp = ncomp_use, intercept = TRUE))
      
      if (nonneg) {
        # Post-hoc non-negativity: clip negative coefficients to zero
        intercept[i] <- b[1]
        coef_values <- b[2:(n_orig + 1)]
        coef_values[coef_values < 0] <- 0
        P[low:high, i] <- coef_values
      } else {
        intercept[i] <- b[1]
        P[low:high, i] <- b[2:(n_orig + 1)]
      }
      
    } else if (model_type %in% c("Ridge", "Lasso")) {
      alpha <- if (model_type == "Ridge") 0 else 1
      
      cvfit <- cv.glmnet(Xmat, yvec,
                         alpha = alpha,
                         intercept = TRUE,
                         standardize = FALSE)
      b <- as.numeric(coef(cvfit, s = "lambda.min"))
      intercept[i] <- b[1]
      
      # Extract only the linear coefficients (first n_orig coefficients after intercept)
      P[low:high, i] <- b[2:(n_orig + 1)]
      
    } else if (model_type == "OLS") {
      X <- cbind(1, Xmat)
      
      # Check for singularity and use appropriate method
      if (ncol(X) >= nrow(X)) {
        # More variables than observations - use ridge penalty
        lambda <- 1e-4
        XtX <- t(X) %*% X
        diag(XtX) <- diag(XtX) + lambda
        b <- as.numeric(solve(XtX, t(X) %*% yvec))
        
      } else {
        # Use QR decomposition (numerically stable)
        qr_decomp <- qr(X)
        
        if (qr_decomp$rank < ncol(X)) {
          # Rank deficient - add small ridge penalty
          lambda <- 1e-6
          XtX <- t(X) %*% X
          diag(XtX) <- diag(XtX) + lambda
          b <- as.numeric(solve(XtX, t(X) %*% yvec))
        } else {
          # Full rank - use standard QR solution
          b <- as.numeric(qr.coef(qr_decomp, yvec))
        }
      }
      
      intercept[i] <- b[1]
      P[low:high, i] <- b[2:(n_orig + 1)]
      
    } else if (model_type == "SVM") {
      if (nonneg) {
        warning("Non-negativity constraint not directly supported for SVM. Proceeding without constraint.")
      }
      
      # CRITICAL: Scale inputs for SVM
      Xmat_scaled <- scale(Xmat)
      yvec_scaled <- scale(yvec)
      
      # Store scaling parameters for back-transformation
      x_center <- attr(Xmat_scaled, "scaled:center")
      x_scale <- attr(Xmat_scaled, "scaled:scale")
      y_center <- attr(yvec_scaled, "scaled:center")
      y_scale <- attr(yvec_scaled, "scaled:scale")
      
      # Handle cases where scale is zero (constant variable)
      x_scale[x_scale == 0] <- 1
      if (y_scale == 0) y_scale <- 1
      
      if (svm_kernel == "linear" || (svm_kernel == "polynomial" && degree == 1)) {
        # Linear SVM - exact coefficients available
        svm_model <- svm(x = Xmat_scaled, y = as.vector(yvec_scaled), 
                         kernel = "linear",
                         type = "eps-regression",
                         cost = 1, 
                         epsilon = 0.1,
                         scale = FALSE)  # Already scaled
        
        # Extract coefficients
        sv_coef <- t(svm_model$coefs) %*% svm_model$SV
        
        # Transform back to original scale
        coef_original <- as.numeric(sv_coef) * (y_scale / x_scale)
        intercept[i] <- y_center - sum(coef_original * x_center) - svm_model$rho * y_scale
        P[low:high, i] <- coef_original
        
      } else if (svm_kernel == "polynomial") {
        # Polynomial SVM with degree > 1
        svm_model <- svm(x = Xmat_scaled, y = as.vector(yvec_scaled), 
                         kernel = "polynomial",
                         type = "eps-regression",
                         degree = degree,
                         coef0 = 1,
                         cost = 1, 
                         gamma = 1 / ncol(Xmat_scaled),
                         epsilon = 0.1,
                         scale = FALSE)
        
        # Approximate linear coefficients using finite differences
        pred_center <- predict(svm_model, Xmat_scaled)
        epsilon_fd <- 1e-5
        approx_coef <- numeric(n_orig)
        
        for (j in 1:n_orig) {
          Xmat_perturb <- Xmat_scaled
          Xmat_perturb[, j] <- Xmat_perturb[, j] + epsilon_fd
          pred_perturb <- predict(svm_model, Xmat_perturb)
          approx_coef[j] <- mean((pred_perturb - pred_center) / epsilon_fd)
        }
        
        # Transform back to original scale
        coef_original <- approx_coef * (y_scale / x_scale)
        intercept[i] <- y_center - sum(coef_original * x_center)
        P[low:high, i] <- coef_original
        
      } else if (svm_kernel == "radial") {
        # Radial basis function (RBF/Gaussian) kernel
        svm_model <- svm(x = Xmat_scaled, y = as.vector(yvec_scaled), 
                         kernel = "radial",
                         type = "eps-regression",
                         cost = 1, 
                         gamma = 1 / ncol(Xmat_scaled),
                         epsilon = 0.1,
                         scale = FALSE)
        
        # Approximate linear coefficients using finite differences
        pred_center <- predict(svm_model, Xmat_scaled)
        epsilon_fd <- 1e-5
        approx_coef <- numeric(n_orig)
        
        for (j in 1:n_orig) {
          Xmat_perturb <- Xmat_scaled
          Xmat_perturb[, j] <- Xmat_perturb[, j] + epsilon_fd
          pred_perturb <- predict(svm_model, Xmat_perturb)
          approx_coef[j] <- mean((pred_perturb - pred_center) / epsilon_fd)
        }
        
        # Transform back to original scale
        coef_original <- approx_coef * (y_scale / x_scale)
        intercept[i] <- y_center - sum(coef_original * x_center)
        P[low:high, i] <- coef_original
        
      } else {
        stop("svm_kernel must be 'linear', 'polynomial', or 'radial'")
      }
      
    } else {
      stop("Unsupported model_type: ", model_type)
    }
    
    pb$tick()
  }
  
  sparse_p <- Matrix(P, sparse = TRUE)
  intcpt <- intercept
  
  assign("sparse_p", sparse_p, envir = .GlobalEnv)
  assign("intcpt", intcpt, envir = .GlobalEnv)
  
  # Apply PDS correction if child.mod exists in global environment
  if (exists("child.mod", envir = .GlobalEnv) && exists("class_val", envir = .GlobalEnv)) {
    child.mod <- get("child.mod", envir = .GlobalEnv)
    class_val <- get("class_val", envir = .GlobalEnv)
    
    child.mod_cor <- as.matrix(child.mod[, -ncol(child.mod), drop = FALSE])
    child.mod_cor <- child.mod_cor %*% as.matrix(sparse_p) + matrix(intcpt, nrow = nrow(child.mod_cor), ncol = length(intcpt), byrow = TRUE)
    child.mod_cor <- as.data.frame(child.mod_cor)
    
    # Determine column names based on whether peaks exists
    if (exists("peaks", envir = .GlobalEnv) && exists("wl.fs", envir = .GlobalEnv)) {
      wl.fs <- get("wl.fs", envir = .GlobalEnv)
      colnames(child.mod_cor) <- wl.fs
    } else if (exists("wl", envir = .GlobalEnv)) {
      wl <- get("wl", envir = .GlobalEnv)
      colnames(child.mod_cor) <- wl
    } else {
      colnames(child.mod_cor) <- colnames(child.mod)[1:(ncol(child.mod) - 1)]
    }
    
    child.mod_cor$Class <- class_val
    
    assign("child.mod_cor", child.mod_cor, envir = .GlobalEnv)
  }
  
  invisible(list(P = sparse_p, intercept = intcpt))
}