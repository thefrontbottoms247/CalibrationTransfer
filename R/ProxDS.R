ProxDS <- function(child_trans, parent_trans, child_mod, lambda = 1e-3) {
  # --- Input prep ---
  child  <- as.matrix(child_trans)
  parent <- as.matrix(parent_trans)
  child_mod <- as.matrix(child_mod)
  n <- nrow(child); m <- ncol(child)
  
  # --- Step 1: Fast distance computation ---
  message("Computing pairwise distances among calibration (child) samples (vectorized)...")
  D <- as.matrix(dist(child, method = "euclidean", upper = TRUE, diag = TRUE))
  
  # --- Step 2: Gaussian weights ---
  message("Computing Gaussian weights...")
  sigma <- mean(D, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) stop("Invalid sigma computed.")
  W_all <- exp(-(D^2) / (2 * sigma^2))
  W_all <- W_all / rowSums(W_all)
  
  # --- Step 3: Parallel correction loop ---
  message("Computing and applying local transfer matrices in parallel...")
  library(doParallel)
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  child_mod_cor <- foreach(i = seq_len(nrow(child_mod)), .combine = rbind, .packages = "base") %dopar% {
    # --- Nearest-neighbor selection ---
    if (nrow(child_mod) == n) {
      idx <- i
    } else {
      dists <- sqrt(rowSums((child - matrix(child_mod[i, ], n, m, byrow = TRUE))^2))
      idx <- which.min(dists)
    }
    
    # --- Local weighting ---
    w <- W_all[idx, ]
    sw <- sqrt(w)
    Xc_w <- child  * sw
    Xp_w <- parent * sw
    
    A <- crossprod(Xc_w) + diag(lambda, m)
    B <- crossprod(Xc_w, Xp_w)
    
    # Numerical sanity checks
    if (any(!is.finite(A)) || any(!is.finite(B))) return(rep(NA_real_, m))
    
    Ti <- tryCatch(solve(A, B),
                   error = function(e) try(qr.solve(A, B), silent = TRUE))
    if (!is.matrix(Ti) || any(!is.finite(Ti))) return(rep(NA_real_, m))
    
    out <- as.numeric(child_mod[i, ]) %*% Ti
    if (any(!is.finite(out))) return(rep(NA_real_, m))
    out
  }
  stopCluster(cl)
  
  # --- Step 4: Results ---
  message("Assigning child.mod_cor to global environment...")
  assign("child.mod_cor", child_mod_cor, envir = .GlobalEnv)
  
  message("ProxDS completed with ",
          sum(!complete.cases(child_mod_cor)),
          " incomplete rows out of ", nrow(child_mod_cor),
          " (parallelized across ", n_cores, " cores).")
  
  invisible(child_mod_cor)
}
