#' @title Nonlinear Piecewise Direct Standardization (N-PDS)
#' @description
#' Extends classical Piecewise Direct Standardization (PDS) by allowing
#' nonlinear, regularized, and machine-learning-based regressions within each
#' local wavelength window. Supports polynomial feature expansion, SVM kernels,
#' and non-negativity constraints.
#'
#' @details
#' When `model_type` is one of "OLS", "Ridge", "Lasso", "PLSR", or "PCR",
#' this function behaves as classical PDS. When using `degree > 1` or
#' `model_type = 'SVM'`, it performs nonlinear piecewise regression approximating
#' the parent–child transfer function.
#'
#' @param parent Numeric matrix or data frame of reference spectra (n × p).
#' @param child  Numeric matrix or data frame of target spectra (n × p).
#' @param child.mod Optional data frame of child spectra to be corrected,
#'        with a `Class` column as the last variable.
#' @param window Character; "static" or "dynamic" window definition method.
#' @param model_type Character; regression type ("OLS", "Ridge", "Lasso",
#'        "PLSR", "PCR", or "SVM").
#' @param degree Integer; polynomial degree for nonlinear expansion (default = 1).
#' @param svm_kernel Character; SVM kernel type ("linear", "polynomial", "radial").
#' @param base_window Integer; base window size for static or non-peak regions.
#' @param peak_window Integer; window size for dynamic peak regions.
#' @param ncomp Integer; number of components for latent variable models.
#' @param ncomp_svd Optional numeric vector of component counts per variable.
#' @param nonneg Logical; whether to enforce non-negativity on regression coefficients.
#'
#' @return A list containing:
#' \item{P}{Sparse coefficient matrix (p × p).}
#' \item{intercept}{Intercept vector (length p).}
#' \item{child.mod_cor}{Corrected dataset if `child.mod` supplied.}
#'
#' @seealso [PDS()] for the strictly linear implementation.
#' @export
NPDS <- function(parent, child, child.mod = NULL,
                 window = "static", model_type = "OLS",
                 degree = 1, svm_kernel = "linear",
                 base_window = 5, peak_window = 20,
                 ncomp = 20, ncomp_svd = NULL, nonneg = FALSE) {
  
  library(progress)
  library(pls)
  library(glmnet)
  library(Matrix)
  library(e1071)
  library(nnls)
  
  stopifnot(ncol(parent) == ncol(child))
  nvar <- ncol(parent)
  
  # --- Determine window sizes ---
  if (window == "dynamic") {
    if (!exists("peaks", envir = .GlobalEnv))
      stop("Dynamic window requires 'peaks' vector in the global environment.")
    peaks <- get("peaks", envir = .GlobalEnv)
    window_n <- rep(base_window, nvar)
    window_n[peaks] <- peak_window
    cat(sprintf("Dynamic window: %d base, %d peak (%d total variables)\n",
                nvar - length(peaks), length(peaks), nvar))
  } else if (window == "static") {
    window_n <- rep(base_window, nvar)
  } else {
    stop("window must be 'static' or 'dynamic'.")
  }
  
  # --- Initialize coefficient and intercept matrices ---
  P <- matrix(0, nrow = nvar, ncol = nvar)
  intercept <- numeric(nvar)
  
  pb <- progress_bar$new(
    format = "  \033[1;32mNPDS building [:bar] :percent :elapsed\033[0m",
    total = nvar, clear = FALSE, width = 70
  )
  
  # --- Main loop ---
  for (i in seq_len(nvar)) {
    half <- floor(window_n[i] / 2)
    low  <- max(1, i - half)
    high <- min(nvar, i + half)
    
    Xmat <- as.matrix(child[, low:high, drop = FALSE])
    yvec <- as.numeric(parent[, i])
    n_orig <- ncol(Xmat)
    
    # Polynomial feature expansion (nonlinear mapping)
    if (degree > 1 && model_type != "SVM") {
      Xmat_poly <- Xmat
      for (d in 2:degree) Xmat_poly <- cbind(Xmat_poly, Xmat^d)
      Xmat <- Xmat_poly
    }
    
    # --- Model fitting ---
    if (nonneg && model_type %in% c("OLS", "Ridge", "Lasso")) {
      if (model_type == "OLS") {
        nnls_fit <- nnls(Xmat, yvec)
        b <- coef(nnls_fit)
        intercept[i] <- mean(yvec - Xmat %*% b)
        P[low:high, i] <- b[1:n_orig]
      } else {
        alpha <- if (model_type == "Ridge") 0 else 1
        cvfit <- cv.glmnet(Xmat, yvec,
                           alpha = alpha,
                           intercept = TRUE,
                           standardize = FALSE,
                           lower.limits = 0)
        b <- as.numeric(coef(cvfit, s = "lambda.min"))
        intercept[i] <- b[1]
        P[low:high, i] <- b[2:(n_orig + 1)]
      }
      
    } else if (model_type %in% c("PLSR", "PCR")) {
      df <- data.frame(y = yvec, Xmat)
      colnames(df) <- c("y", paste0("X", seq_len(ncol(Xmat))))
      ncomp_use <- if (!is.null(ncomp_svd))
        min(ncomp_svd[i], ncol(Xmat)) else min(ncomp, ncol(Xmat))
      
      fit <- if (model_type == "PLSR") {
        plsr(y ~ ., data = df, ncomp = ncomp_use, validation = "none")
      } else {
        pcr(y ~ ., data = df, ncomp = ncomp_use, validation = "none")
      }
      
      b <- drop(coef(fit, ncomp = ncomp_use, intercept = TRUE))
      intercept[i] <- b[1]
      coef_values <- b[2:(n_orig + 1)]
      if (nonneg) coef_values[coef_values < 0] <- 0
      P[low:high, i] <- coef_values
      
    } else if (model_type %in% c("Ridge", "Lasso")) {
      alpha <- if (model_type == "Ridge") 0 else 1
      cvfit <- cv.glmnet(Xmat, yvec,
                         alpha = alpha,
                         intercept = TRUE,
                         standardize = FALSE)
      b <- as.numeric(coef(cvfit, s = "lambda.min"))
      intercept[i] <- b[1]
      P[low:high, i] <- b[2:(n_orig + 1)]
      
    } else if (model_type == "OLS") {
      X <- cbind(1, Xmat)
      XtX <- crossprod(X)
      if (ncol(X) >= nrow(X) || qr(X)$rank < ncol(X)) {
        lambda <- 1e-4
        diag(XtX) <- diag(XtX) + lambda
        b <- as.numeric(solve(XtX, crossprod(X, yvec)))
      } else {
        b <- as.numeric(qr.coef(qr(X), yvec))
      }
      intercept[i] <- b[1]
      P[low:high, i] <- b[2:(n_orig + 1)]
      
    } else if (model_type == "SVM") {
      if (nonneg) warning("Non-negativity not supported for SVM. Proceeding unconstrained.")
      Xmat_scaled <- scale(Xmat)
      yvec_scaled <- scale(yvec)
      x_center <- attr(Xmat_scaled, "scaled:center")
      x_scale  <- attr(Xmat_scaled, "scaled:scale")
      y_center <- attr(yvec_scaled, "scaled:center")
      y_scale  <- attr(yvec_scaled, "scaled:scale")
      x_scale[x_scale == 0] <- 1
      if (y_scale == 0) y_scale <- 1
      
      if (svm_kernel == "linear" || (svm_kernel == "polynomial" && degree == 1)) {
        svm_model <- svm(x = Xmat_scaled, y = yvec_scaled,
                         kernel = "linear", type = "eps-regression",
                         cost = 1, epsilon = 0.1, scale = FALSE)
        sv_coef <- t(svm_model$coefs) %*% svm_model$SV
        coef_original <- as.numeric(sv_coef) * (y_scale / x_scale)
        intercept[i] <- y_center - sum(coef_original * x_center) - svm_model$rho * y_scale
        P[low:high, i] <- coef_original
      } else {
        # Nonlinear kernels: finite-difference approximation
        svm_model <- svm(x = Xmat_scaled, y = yvec_scaled,
                         kernel = svm_kernel, type = "eps-regression",
                         degree = degree, coef0 = 1,
                         cost = 1, gamma = 1 / ncol(Xmat_scaled),
                         epsilon = 0.1, scale = FALSE)
        pred_center <- predict(svm_model, Xmat_scaled)
        epsilon_fd <- 1e-5
        approx_coef <- numeric(n_orig)
        for (j in seq_len(n_orig)) {
          Xp <- Xmat_scaled
          Xp[, j] <- Xp[, j] + epsilon_fd
          pred_p <- predict(svm_model, Xp)
          approx_coef[j] <- mean((pred_p - pred_center) / epsilon_fd)
        }
        coef_original <- approx_coef * (y_scale / x_scale)
        intercept[i] <- y_center - sum(coef_original * x_center)
        P[low:high, i] <- coef_original
      }
      
    } else {
      stop("Unsupported model_type: ", model_type)
    }
    
    pb$tick()
  }
  
  # --- Package results ---
  sparse_p <- Matrix(P, sparse = TRUE)
  intcpt <- intercept
  
  assign("sparse_p", sparse_p, envir = .GlobalEnv)
  assign("intcpt", intcpt, envir = .GlobalEnv)
  
  # --- Apply correction if child.mod provided ---
  if (!is.null(child.mod)) {
    Xmod <- as.matrix(child.mod[, -ncol(child.mod), drop = FALSE])
    child.mod_cor <- Xmod %*% as.matrix(sparse_p) +
      matrix(intcpt, nrow = nrow(Xmod), ncol = length(intcpt), byrow = TRUE)
    child.mod_cor <- as.data.frame(child.mod_cor)
    
    # Column names from parent if available
    if (!is.null(colnames(parent))) {
      colnames(child.mod_cor) <- colnames(parent)
    } else {
      colnames(child.mod_cor) <- colnames(child.mod)[1:(ncol(child.mod) - 1)]
    }
    
    # Preserve Class column if present
    if ("Class" %in% colnames(child.mod))
      child.mod_cor$Class <- child.mod$Class
    
    assign("child.mod_cor", child.mod_cor, envir = .GlobalEnv)
  }
  
  invisible(list(P = sparse_p, intercept = intcpt))
}
