#' @title Stabilized Thin-Plate Dynamic PCA Correction
#' @description
#' Corrects child PCA scores to align with parent PCA using hyperspherical
#' coordinates and thin-plate spline smoothing. Supports arbitrary dimensions.
#'
#' @param pca_parent Parent PCA object (prcomp)
#' @param pca_child Child PCA object (prcomp)
#' @param ncomp Number of principal components to correct (default = 2)
#' @param lambda Smoothing parameter for thin-plate splines (default = 0.01)
#' @param radial_quantile Quantile for radial damping threshold (default = 0.85)
#' @param angle_quantile Quantile for angular damping threshold (default = 0.765)
#' @param pc1_anisotropy Strength of PC1 anisotropic damping (default = 0.5)
#' @param blend_scale Scaling factor for adaptive blending (default = 1.5)
#' @param max_blend Maximum blending weight (default = 0.4)
#' @param min_blend Minimum blending weight (default = 0.05)
#'
#' @return List containing:
#'   \item{child_corrected}{Corrected child PCA scores}
#'   \item{diagnostics}{Data frame with intermediate values for diagnostics}
#'   \item{models}{List of fitted TPS models for each coordinate}
#'
#' @examples
#' \dontrun{
#' library(fields)
#' result <- pca_tps_correct(parent.pca, child.pca, ncomp = 3)
#' child.pca_cor <- result$child_corrected
#' }
#'
#' @export
PCA_Correction <- function(pca_parent, 
                            pca_child, 
                            ncomp = 2,
                            lambda = 0.01,
                            radial_quantile = 0.85,
                            angle_quantile = 0.765,
                            pc1_anisotropy = 0.5,
                            blend_scale = 1.5,
                            max_blend = 0.4,
                            min_blend = 0.05) {
  
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Package 'fields' is required. Install it with: install.packages('fields')")
  }
  
  # --- Extract scores ---
  scores_parent <- as.matrix(pca_parent$x[, 1:ncomp])
  scores_child <- as.matrix(pca_child$x[, 1:ncomp])
  
  if (nrow(scores_parent) != nrow(scores_child)) {
    stop("Parent and child must have the same number of samples.")
  }
  
  n_samples <- nrow(scores_child)
  n_angles <- ncomp - 1
  
  # --- Centroids ---
  centroid_parent <- colMeans(scores_parent)
  centroid_child <- colMeans(scores_child)
  
  # --- Center scores ---
  centered_parent <- sweep(scores_parent, 2, centroid_parent, "-")
  centered_child <- sweep(scores_child, 2, centroid_child, "-")
  
  # ============================================================
  # Convert to hyperspherical coordinates
  # ============================================================
  
  # Radial distance
  r_parent <- sqrt(rowSums(centered_parent^2))
  r_child <- sqrt(rowSums(centered_child^2))
  
  # Angular coordinates
  angles_parent <- matrix(NA, nrow = n_samples, ncol = n_angles)
  angles_child <- matrix(NA, nrow = n_samples, ncol = n_angles)
  
  for (i in seq_len(n_samples)) {
    # Parent angles
    u_p <- centered_parent[i, ]
    r_p <- r_parent[i]
    if (r_p > 1e-10) {
      u_p <- u_p / r_p
      for (k in seq_len(n_angles - 1)) {
        denom <- sqrt(sum(u_p[k:ncomp]^2))
        angles_parent[i, k] <- acos(pmax(pmin(u_p[k] / denom, 1), -1))
      }
      angles_parent[i, n_angles] <- atan2(u_p[ncomp], u_p[ncomp - 1])
    }
    
    # Child angles
    u_c <- centered_child[i, ]
    r_c <- r_child[i]
    if (r_c > 1e-10) {
      u_c <- u_c / r_c
      for (k in seq_len(n_angles - 1)) {
        denom <- sqrt(sum(u_c[k:ncomp]^2))
        angles_child[i, k] <- acos(pmax(pmin(u_c[k] / denom, 1), -1))
      }
      angles_child[i, n_angles] <- atan2(u_c[ncomp], u_c[ncomp - 1])
    }
  }
  
  # ============================================================
  # Prepare predictors and targets
  # ============================================================
  
  # Create predictor matrix: [r_child, sin/cos of all child angles]
  predictors <- cbind(r_child)
  for (k in seq_len(n_angles)) {
    predictors <- cbind(predictors, sin(angles_child[, k]), cos(angles_child[, k]))
  }
  
  # Target: radial ratio
  r_ratio <- r_parent / pmax(r_child, 1e-6)
  
  # Targets: angular differences (wrapped)
  angle_diffs <- matrix(NA, nrow = n_samples, ncol = n_angles)
  for (k in seq_len(n_angles)) {
    angle_diffs[, k] <- atan2(
      sin(angles_parent[, k] - angles_child[, k]),
      cos(angles_parent[, k] - angles_child[, k])
    )
  }
  
  # ============================================================
  # Fit thin-plate splines
  # ============================================================
  
  message("Fitting thin-plate spline models...")
  
  model_ratio <- fields::Tps(predictors, r_ratio, lambda = lambda)
  pred_ratio <- predict(model_ratio, predictors)
  
  models_angle <- vector("list", n_angles)
  pred_angles <- matrix(NA, nrow = n_samples, ncol = n_angles)
  
  for (k in seq_len(n_angles)) {
    models_angle[[k]] <- fields::Tps(predictors, angle_diffs[, k], lambda = lambda)
    pred_angles[, k] <- predict(models_angle[[k]], predictors)
  }
  
  # ============================================================
  # Controlled corrections with damping
  # ============================================================
  
  q_radial <- quantile(r_child, radial_quantile)
  q_angle <- quantile(r_child, angle_quantile)
  sd_r <- sd(r_child)
  
  # PC1 weight for anisotropic damping (stronger damping along PC1)
  pc1_weight <- abs(cos(angles_child[, 1]))
  
  # Radial damping
  radial_damp <- 1 / (1 + exp((r_child - q_radial) / sd_r)) * 
    (1 - pc1_anisotropy * pc1_weight)
  
  # Angular damping
  angle_damp <- 1 / (1 + exp((r_child - q_angle) / (0.5 * sd_r)))
  
  # Apply damped corrections
  r_hat <- r_child * (1 + (pred_ratio - 1) * radial_damp)
  angles_hat <- angles_child + pred_angles * angle_damp
  
  # ============================================================
  # Adaptive blending based on prediction error
  # ============================================================
  
  err_r <- abs(r_ratio - pred_ratio)
  loc_blend <- pmin(max_blend, min_blend + blend_scale * err_r / max(err_r, na.rm = TRUE))
  
  r_hat <- (1 - loc_blend) * r_hat + loc_blend * r_child
  angles_hat <- (1 - loc_blend) * angles_hat + loc_blend * angles_child
  
  # ============================================================
  # Convert back to Cartesian coordinates
  # ============================================================
  
  child_corrected <- matrix(NA, nrow = n_samples, ncol = ncomp)
  
  for (i in seq_len(n_samples)) {
    r <- r_hat[i]
    angles <- angles_hat[i, ]
    
    # Reconstruct Cartesian from hyperspherical
    coords <- numeric(ncomp)
    sin_prod <- 1
    
    for (k in seq_len(n_angles - 1)) {
      coords[k] <- sin_prod * cos(angles[k])
      sin_prod <- sin_prod * sin(angles[k])
    }
    
    coords[ncomp - 1] <- sin_prod * cos(angles[n_angles])
    coords[ncomp] <- sin_prod * sin(angles[n_angles])
    
    child_corrected[i, ] <- r * coords
  }
  
  # Add back parent centroid
  child_corrected <- sweep(child_corrected, 2, centroid_parent, "+")
  
  # Final recentering
  offset <- colMeans(child_corrected) - centroid_parent
  child_corrected <- sweep(child_corrected, 2, offset, "-")
  
  # ============================================================
  # Prepare output
  # ============================================================
  
  colnames(child_corrected) <- colnames(scores_child)
  rownames(child_corrected) <- rownames(scores_child)
  
  diagnostics <- data.frame(
    r_child = r_child,
    r_parent = r_parent,
    r_hat = r_hat,
    r_ratio = r_ratio,
    pred_ratio = pred_ratio,
    radial_damp = radial_damp,
    angle_damp = angle_damp,
    blend_weight = loc_blend
  )
  
  message("✅ PCA correction complete!")
  
  return(list(
    child_corrected = child_corrected,
    diagnostics = diagnostics,
    models = list(
      ratio = model_ratio,
      angles = models_angle
    )
  ))
}


#' @title Visualize PCA Correction Results
#' @description
#' Creates diagnostic plots for PCA correction in 2D or 3D.
#'
#' @param pca_parent Parent PCA object
#' @param pca_child Child PCA object
#' @param result Output from pca_tps_correct
#' @param dims Dimensions to plot (default = c(1,2))
#'
#' @export
plot_pca_correction <- function(pca_parent, pca_child, result, dims = c(1, 2)) {
  
  if (length(dims) > 2) {
    warning("Only 2D plotting supported. Using first two dimensions specified.")
    dims <- dims[1:2]
  }
  
  pc1 <- dims[1]
  pc2 <- dims[2]
  
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  
  plot(pca_parent$x[, pc1], pca_parent$x[, pc2], 
       pch = 16, col = "firebrick", cex = 0.8,
       xlab = paste0("PC", pc1), 
       ylab = paste0("PC", pc2),
       main = "Thin-Plate Spline PCA Correction")
  
  points(pca_child$x[, pc1], pca_child$x[, pc2], 
         pch = 16, col = "darkgreen", cex = 0.8)
  
  points(result$child_corrected[, pc1], result$child_corrected[, pc2], 
         pch = 16, col = "navy", cex = 0.8)
  
  legend("topright",
         legend = c("Parent", "Child", "Child Corrected"),
         col = c("firebrick", "darkgreen", "navy"), 
         pch = 16, cex = 0.9)
}