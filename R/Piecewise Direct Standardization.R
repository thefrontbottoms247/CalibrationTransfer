#' @title Piecewise Direct Standardization (PDS)
#' @description
#' Performs classical piecewise direct standardization between parent and child
#' instruments. Each wavelength region in the child is locally regressed onto
#' the corresponding region in the parent using linear models.
#'
#' @param parent Numeric matrix or data frame of reference spectra.
#' @param child  Numeric matrix or data frame of target spectra.
#' @param child.mod Optional data frame of child spectra to be corrected,
#'        with a `Class` column as the last variable.
#' @param window Character; "static" or "dynamic" window definition method.
#' @param ncomp Integer; number of components for latent variable models.
#' @param model_type Character; regression type ("OLS", "Ridge", "Lasso",
#'        "PLSR", or "PCR").
#' @param base_window Integer; base window size for static or non-peak regions.
#' @param peak_window Integer; window size for dynamic peak regions.
#' @param ncomp_svd Optional numeric vector of component counts per variable.
#' @param nonneg Logical; whether to enforce non-negativity on regression coefficients.
#'
#' @return A list containing:
#'   \item{P}{Sparse coefficient matrix (p × p).}
#'   \item{intercept}{Intercept vector (length p).}
#'   \item{child.mod_cor}{Corrected dataset (if `child.mod` supplied).}
#' @export
PDS <- function(parent, child, child.mod,
                window = "static", ncomp = 20,
                model_type = "OLS", base_window = 5,
                peak_window = 20, ncomp_svd = NULL,
                nonneg = FALSE) {

  library(progress)
  library(pls)
  library(glmnet)
  library(Matrix)
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
    cat(sprintf("Dynamic window: %d base, %d peak (%d total)\n",
                nvar - length(peaks), length(peaks), nvar))
  } else if (window == "static") {
    window_n <- rep(base_window, nvar)
  } else {
    stop("window must be 'static' or 'dynamic'.")
  }

  # --- Initialize ---
  P <- matrix(0, nrow = nvar, ncol = nvar)
  intercept <- numeric(nvar)

  pb <- progress_bar$new(
    format = "  \033[1;32mPDS building [:bar] :percent :elapsed\033[0m",
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

    # Determine number of components
    ncomp_use <- ncomp
    if (!is.null(ncomp_svd) && model_type %in% c("PLSR", "PCR")) {
      if (length(ncomp_svd) != nvar)
        stop("ncomp_svd length must equal number of variables (", nvar, ")")
      ncomp_use <- min(ncomp_svd[i], ncol(Xmat))
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
        cvfit <- cv.glmnet(Xmat, yvec, alpha = alpha,
                           intercept = TRUE, standardize = FALSE,
                           lower.limits = 0)
        b <- as.numeric(coef(cvfit, s = "lambda.min"))
        intercept[i] <- b[1]
        P[low:high, i] <- b[2:(n_orig + 1)]
      }

    } else if (model_type %in% c("PLSR", "PCR")) {
      df <- data.frame(y = yvec, Xmat)
      colnames(df) <- c("y", paste0("X", seq_len(ncol(Xmat))))
      fit <- if (model_type == "PLSR") {
        plsr(y ~ ., data = df, ncomp = ncomp_use, validation = "none")
      } else {
        pcr(y ~ ., data = df, ncomp = ncomp_use, validation = "none")
      }
      b <- drop(coef(fit, ncomp = ncomp_use, intercept = TRUE))
      intercept[i] <- b[1]
      coefs <- b[2:(n_orig + 1)]
      if (nonneg) coefs[coefs < 0] <- 0
      P[low:high, i] <- coefs

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

    } else {
      stop("Unsupported model_type: ", model_type)
    }

    pb$tick()
  }

  # --- Finalize ---
  sparse_p <- Matrix(P, sparse = TRUE)
  intcpt <- intercept
  assign("sparse_p", sparse_p, envir = .GlobalEnv)
  assign("intcpt", intcpt, envir = .GlobalEnv)

  # --- Apply to child.mod if provided ---
  if (!is.null(child.mod)) {
    Xmod <- as.matrix(child.mod[, -ncol(child.mod), drop = FALSE])
    child.mod_cor <- Xmod %*% as.matrix(sparse_p) +
      matrix(intcpt, nrow = nrow(Xmod), ncol = length(intcpt), byrow = TRUE)
    child.mod_cor <- as.data.frame(child.mod_cor)

    if (!is.null(colnames(parent))) {
      colnames(child.mod_cor) <- colnames(parent)
    } else {
      colnames(child.mod_cor) <- colnames(child.mod)[1:(ncol(child.mod) - 1)]
    }

    if ("Class" %in% colnames(child.mod))
      child.mod_cor$Class <- child.mod$Class

    assign("child.mod_cor", child.mod_cor, envir = .GlobalEnv)
  }

  invisible(list(P = sparse_p, intercept = intcpt))
}
