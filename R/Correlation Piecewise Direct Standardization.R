#' @title CPDS: Combined (C)ORAL and (P)iecewise (D)irect (S)tandardization
#' @description
#' Flexible pipeline for calibration transfer:
#' - order = "PDS->CORAL": run PDS, then global CORAL on the PDS output.
#' - order = "CORAL->PDS": global CORAL first, then PDS.
#' - order = "local-CORAL-PDS": CORAL whitening/recoloring inside each PDS window.
#'
#' @param parent matrix/data.frame; reference spectra.
#' @param child  matrix/data.frame; target spectra to learn mapping from.
#' @param child.mod optional data frame of child spectra to correct; last column may be Class.
#' @param order one of c("PDS->CORAL","CORAL->PDS","local-CORAL-PDS").
#' @param window "static" or "dynamic" PDS windows. For "dynamic" a global `peaks` integer vector must exist.
#' @param base_window integer; base PDS window size.
#' @param peak_window integer; dynamic peak window size.
#' @param model_type "OLS","Ridge","Lasso","PLSR","PCR".
#' @param ncomp integer; components for PLSR/PCR.
#' @param nonneg logical; enforce non-negativity in linear models where applicable.
#' @param coral_method "auto","full","svd".
#' @param coral_center logical; CORAL centering.
#' @param coral_scale logical; CORAL variance scaling.
#' @param coral_lambda numeric; CORAL ridge.
#' @return list with elements:
#'   \item{X_child_aligned}{transformed child (training) spectra}
#'   \item{P}{sparse PDS coefficient matrix (if PDS used)}
#'   \item{intercept}{PDS intercept vector (if PDS used)}
#'   \item{A}{global CORAL transform (if used)}
#'   \item{child.mod_cor}{corrected child.mod if supplied}
#' @export
CPDS <- function(parent, child, child.mod = NULL,
                 order = c("PDS->CORAL","CORAL->PDS","local-CORAL-PDS"),
                 window = "static", base_window = 5, peak_window = 20,
                 model_type = "OLS", ncomp = 20, nonneg = FALSE,
                 coral_method = "auto", coral_center = TRUE, coral_scale = FALSE,
                 coral_lambda = 1e-6) {
  
  order <- match.arg(order)
  parent <- as.matrix(parent); child <- as.matrix(child)
  stopifnot(ncol(parent) == ncol(child))
  p <- ncol(parent)
  
  suppressPackageStartupMessages({
    require(Matrix); require(glmnet); require(pls); require(nnls); require(MASS); require(expm)
  })
  
  .coral_fit <- function(Xc, Xp, center = TRUE, scale = FALSE, lambda = 1e-6, method = c("auto","full","svd")) {
    method <- match.arg(method)
    Xc <- as.matrix(Xc); Xp <- as.matrix(Xp)
    nc <- nrow(Xc); np <- nrow(Xp); pp <- ncol(Xc)
    if (ncol(Xp) != pp) stop("CORAL: mismatched columns.")
    if (method == "auto") method <- if (pp > (nc + np)) "svd" else "full"
    
    if (center) {
      muc <- colMeans(Xc); mup <- colMeans(Xp)
      Xc <- sweep(Xc, 2, muc); Xp <- sweep(Xp, 2, mup)
    }
    if (scale) {
      sdc <- pmax(apply(Xc, 2, sd), 1e-12)
      sdp <- pmax(apply(Xp, 2, sd), 1e-12)
      Xc <- sweep(Xc, 2, sdc, "/"); Xp <- sweep(Xp, 2, sdp, "/")
    }
    
    if (method == "full") {
      Sc <- cov(Xc) + diag(lambda, pp)
      Sp <- cov(Xp) + diag(lambda, pp)
      Ac <- MASS::ginv(expm::sqrtm(Sc))
      Ap <- expm::sqrtm(Sp)
      A  <- Ac %*% Ap
    } else {
      svc <- svd(Xc, nu = 0); svp <- svd(Xp, nu = 0)
      r   <- min(length(svc$d), length(svp$d), nc - 1, np - 1)
      Vc  <- svc$v[, 1:r, drop = FALSE]; Vp <- svp$v[, 1:r, drop = FALSE]
      Scd <- svc$d[1:r]; Spd <- svp$d[1:r]
      A   <- Vc %*% (diag(1 / sqrt(pmax(Scd^2, lambda)), r) %*% (t(Vc) %*% Vp) %*%
                       diag(sqrt(pmax(Spd^2, lambda)), r)) %*% t(Vp)
    }
    list(A = A, center = if (center) mup else NULL)
  }
  
  .coral_apply <- function(X, fit) {
    Y <- X %*% fit$A
    if (!is.null(fit$center)) Y <- sweep(Y, 2, fit$center, "+")
    Y
  }
  
  .pds_fit <- function(Parent, Child,
                       window = "static", base_window = 5, peak_window = 20,
                       model_type = "OLS", ncomp = 20, nonneg = FALSE,
                       local_coral = FALSE,
                       coral_lambda = 1e-6) {
    
    nvar <- ncol(Parent)
    if (window == "dynamic") {
      if (!exists("peaks", envir = .GlobalEnv)) stop("Dynamic window needs global `peaks`.")
      pk <- get("peaks", envir = .GlobalEnv)
      wv <- rep(base_window, nvar); wv[pk] <- peak_window
    } else if (window == "static") {
      wv <- rep(base_window, nvar)
    } else stop("window must be 'static' or 'dynamic'.")
    
    P <- matrix(0, nvar, nvar)
    icpt <- numeric(nvar)
    
    for (i in seq_len(nvar)) {
      half <- floor(wv[i] / 2)
      lo <- max(1, i - half); hi <- min(nvar, i + half)
      Xw <- as.matrix(Child[, lo:hi, drop = FALSE])
      yw <- as.numeric(Parent[, i])
      
      if (local_coral && (hi - lo + 1) >= 2) {
        Sc <- cov(Child[, lo:hi, drop = FALSE]) + diag(coral_lambda, (hi - lo + 1))
        Sp <- cov(Parent[, lo:hi, drop = FALSE]) + diag(coral_lambda, (hi - lo + 1))
        A  <- MASS::ginv(expm::sqrtm(Sc)) %*% expm::sqrtm(Sp)
        Xw <- Xw %*% A
      }
      
      if (nonneg && model_type %in% c("OLS","Ridge","Lasso")) {
        if (model_type == "OLS") {
          fit <- nnls(Xw, yw)
          b <- coef(fit); icpt[i] <- mean(yw - Xw %*% b); P[lo:hi, i] <- b
        } else {
          alpha <- if (model_type == "Ridge") 0 else 1
          cf <- cv.glmnet(Xw, yw, alpha = alpha, intercept = TRUE, standardize = FALSE, lower.limits = 0)
          bb <- as.numeric(coef(cf, s = "lambda.min")); icpt[i] <- bb[1]; P[lo:hi, i] <- bb[-1]
        }
      } else if (model_type %in% c("PLSR","PCR")) {
        df <- data.frame(y = yw, Xw); colnames(df) <- c("y", paste0("x", seq_len(ncol(Xw))))
        fit <- if (model_type == "PLSR") pls::plsr(y ~ ., ncomp = min(ncomp, ncol(Xw)), data = df, validation = "none")
        else                      pls::pcr (y ~ ., ncomp = min(ncomp, ncol(Xw)), data = df, validation = "none")
        bb <- drop(coef(fit, ncomp = min(ncomp, ncol(Xw)), intercept = TRUE))
        icpt[i] <- bb[1]; co <- bb[-1]; if (nonneg) co[co < 0] <- 0; P[lo:hi, i] <- co
      } else if (model_type %in% c("Ridge","Lasso")) {
        alpha <- if (model_type == "Ridge") 0 else 1
        cf <- cv.glmnet(Xw, yw, alpha = alpha, intercept = TRUE, standardize = FALSE)
        bb <- as.numeric(coef(cf, s = "lambda.min")); icpt[i] <- bb[1]; P[lo:hi, i] <- bb[-1]
      } else if (model_type == "OLS") {
        X <- cbind(1, Xw); XtX <- crossprod(X)
        if (ncol(X) >= nrow(X) || qr(X)$rank < ncol(X)) diag(XtX) <- diag(XtX) + 1e-4
        b <- as.numeric(solve(XtX, crossprod(X, yw)))
        icpt[i] <- b[1]; P[lo:hi, i] <- b[-1]
      } else stop("Unsupported model_type.")
    }
    
    list(P = Matrix(P, sparse = TRUE), intercept = icpt)
  }
  
  apply_pds <- function(X, fit) {
    X %*% as.matrix(fit$P) + matrix(fit$intercept, nrow(X), p, byrow = TRUE)
  }
  
  # ---- Pipelines ----
  if (order == "CORAL->PDS") {
    cf  <- .coral_fit(child, parent, center = coral_center, scale = coral_scale,
                      lambda = coral_lambda, method = coral_method)
    child_c <- .coral_apply(child, cf)
    pfit <- .pds_fit(parent, child_c, window, base_window, peak_window,
                     model_type, ncomp, nonneg, local_coral = FALSE, coral_lambda = coral_lambda)
    Xc_aligned <- child_c
    A <- cf$A
    
  } else if (order == "PDS->CORAL") {
    pfit <- .pds_fit(parent, child, window, base_window, peak_window,
                     model_type, ncomp, nonneg, local_coral = FALSE, coral_lambda = coral_lambda)
    Xp_child <- apply_pds(child, pfit)
    cf  <- .coral_fit(Xp_child, parent, center = coral_center, scale = coral_scale,
                      lambda = coral_lambda, method = coral_method)
    Xc_aligned <- .coral_apply(Xp_child, cf)
    A <- cf$A
    
  } else { # local-CORAL-PDS
    pfit <- .pds_fit(parent, child, window, base_window, peak_window,
                     model_type, ncomp, nonneg, local_coral = TRUE, coral_lambda = coral_lambda)
    Xc_aligned <- child  # training child in original space; mapping lives in P
    A <- NULL
  }
  
  # ---- child.mod correction ----
  child.mod_cor <- NULL
  if (!is.null(child.mod)) {
    Xm <- as.matrix(child.mod[, -ncol(child.mod), drop = FALSE])
    if (order == "CORAL->PDS") {
      Xm <- .coral_apply(Xm, .coral_fit(child, parent, center = coral_center, scale = coral_scale,
                                        lambda = coral_lambda, method = coral_method))
      child.mod_cor <- apply_pds(Xm, pfit)
    } else if (order == "PDS->CORAL") {
      Xm <- apply_pds(Xm, pfit)
      child.mod_cor <- .coral_apply(Xm, .coral_fit(apply_pds(child, pfit), parent,
                                                   center = coral_center, scale = coral_scale,
                                                   lambda = coral_lambda, method = coral_method))
    } else { # local-CORAL-PDS uses same local windows as training
      child.mod_cor <- apply_pds(as.matrix(child.mod[, -ncol(child.mod), drop = FALSE]), pfit)
    }
    child.mod_cor <- as.data.frame(child.mod_cor)
    if (!is.null(colnames(parent))) colnames(child.mod_cor) <- colnames(parent)
    if ("Class" %in% colnames(child.mod)) child.mod_cor$Class <- child.mod$Class
  }
  
  out <- list(X_child_aligned = Xc_aligned,
              P = if (exists("pfit")) pfit$P else NULL,
              intercept = if (exists("pfit")) pfit$intercept else NULL,
              A = A,
              child.mod_cor = child.mod_cor)
  invisible(out)
}
