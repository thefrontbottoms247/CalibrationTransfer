#' @title Stabilized Tensor-Product PCA Correction
#' @description
#' Corrects child PCA scores to align with parent PCA using hyperspherical
#' coordinates and smooth additive splines (`mgcv::gam`). Handles
#' arbitrary PCA dimensionalities and avoids memory limits.
#' Automatically plots the corrected scores.
#'
#' @param pca_parent Parent PCA object (from prcomp)
#' @param pca_child  Child PCA object (from prcomp)
#' @param ncomp Number of principal components to correct (default = 5)
#' @param lambda Regularization factor for spline smoothing (default = 0.01)
#' @param radial_quantile Quantile for radial damping (default = 0.85)
#' @param angle_quantile Quantile for angular damping (default = 0.765)
#' @param pc1_anisotropy Damping strength along PC1 (default = 0.5)
#' @param blend_scale Scaling for adaptive blending (default = 1.5)
#' @param max_blend Maximum blending weight (default = 0.4)
#' @param min_blend Minimum blending weight (default = 0.05)
#' @param plot_dims PCs to visualize (default = c(1,2))
#' @param plot Whether to produce a diagnostic plot (default = TRUE)
#'
#' @return List with:
#' \itemize{
#'   \item{child_corrected}{Corrected PCA scores for the child set.}
#'   \item{diagnostics}{Data frame of intermediate diagnostics.}
#'   \item{models}{List of fitted GAM models.}
#' }
#'
#' @export
PCA_Correction <- function(pca_parent,
                           pca_child,
                           ncomp = 5,
                           lambda = 0.01,
                           radial_quantile = 0.85,
                           angle_quantile = 0.765,
                           pc1_anisotropy = 0.5,
                           blend_scale = 1.5,
                           max_blend = 0.4,
                           min_blend = 0.05,
                           plot_dims = c(1, 2),
                           plot = TRUE) {

  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Package 'mgcv' is required. Install with install.packages('mgcv').")

  # --- Extract PCA scores ---
  scores_parent <- as.matrix(pca_parent$x[, 1:ncomp])
  scores_child  <- as.matrix(pca_child$x[, 1:ncomp])

  if (nrow(scores_parent) != nrow(scores_child))
    stop("Parent and child must have the same number of rows (samples).")

  n_samples <- nrow(scores_child)
  n_angles  <- ncomp - 1

  # --- Compute centroids ---
  centroid_parent <- colMeans(scores_parent)
  centroid_child  <- colMeans(scores_child)

  # --- Center the PCA scores ---
  centered_parent <- sweep(scores_parent, 2, centroid_parent, "-")
  centered_child  <- sweep(scores_child, 2, centroid_child, "-")

  # ============================================================
  # Convert to hyperspherical coordinates
  # ============================================================

  r_parent <- sqrt(rowSums(centered_parent^2))
  r_child  <- sqrt(rowSums(centered_child^2))

  angles_parent <- matrix(NA, n_samples, n_angles)
  angles_child  <- matrix(NA, n_samples, n_angles)

  for (i in seq_len(n_samples)) {
    # Parent
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

    # Child
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

  # Predictor matrix: radius + sin/cos terms for each angle
  predictors <- cbind(r_child)
  for (k in seq_len(n_angles)) {
    predictors <- cbind(predictors, sin(angles_child[, k]), cos(angles_child[, k]))
  }
  predictors <- scale(predictors)
  pred_df <- as.data.frame(predictors)
  colnames(pred_df) <- paste0("V", seq_len(ncol(pred_df)))

  # Targets
  r_ratio <- r_parent / pmax(r_child, 1e-6)
  angle_diffs <- matrix(NA, n_samples, n_angles)
  for (k in seq_len(n_angles)) {
    angle_diffs[, k] <- atan2(
      sin(angles_parent[, k] - angles_child[, k]),
      cos(angles_parent[, k] - angles_child[, k])
    )
  }

  # ============================================================
  # Fit additive spline models (memory-efficient alternative)
  # ============================================================

  library(mgcv)
  message("Fitting additive spline models...")

  # Build additive formula with smooth terms for each predictor
  smooth_terms <- paste0("s(", colnames(pred_df), ", k=5)", collapse = " + ")
  formula_ratio <- as.formula(paste("r_ratio ~", smooth_terms))

  model_ratio <- mgcv::gam(formula_ratio,
                           data = cbind(pred_df, r_ratio = r_ratio),
                           method = "REML")

  pred_ratio <- predict(model_ratio, newdata = pred_df)

  models_angle <- vector("list", n_angles)
  pred_angles  <- matrix(NA, nrow = n_samples, ncol = n_angles)

  for (k in seq_len(n_angles)) {
    formula_angle <- as.formula(
      paste("angle_diffs[,", k, "] ~", smooth_terms)
    )
    models_angle[[k]] <- mgcv::gam(formula_angle,
                                   data = cbind(pred_df, angle = angle_diffs[, k]),
                                   method = "REML")
    pred_angles[, k] <- predict(models_angle[[k]], newdata = pred_df)
  }

  # ============================================================
  # Controlled damping and blending
  # ============================================================

  q_radial <- quantile(r_child, radial_quantile)
  q_angle  <- quantile(r_child, angle_quantile)
  sd_r <- sd(r_child)

  pc1_weight <- abs(cos(angles_child[, 1]))
  radial_damp <- 1 / (1 + exp((r_child - q_radial) / sd_r)) *
    (1 - pc1_anisotropy * pc1_weight)
  angle_damp <- 1 / (1 + exp((r_child - q_angle) / (0.5 * sd_r)))

  r_hat <- r_child * (1 + (pred_ratio - 1) * radial_damp)
  angles_hat <- angles_child + pred_angles * angle_damp

  err_r <- abs(r_ratio - pred_ratio)
  loc_blend <- pmin(max_blend,
                    min_blend + blend_scale * err_r / max(err_r, na.rm = TRUE))

  r_hat <- (1 - loc_blend) * r_hat + loc_blend * r_child
  angles_hat <- (1 - loc_blend) * angles_hat + loc_blend * angles_child

  # ============================================================
  # Convert back to Cartesian coordinates
  # ============================================================

  child_corrected <- matrix(NA, n_samples, ncomp)

  for (i in seq_len(n_samples)) {
    r <- r_hat[i]
    ang <- angles_hat[i, ]
    coords <- numeric(ncomp)
    sin_prod <- 1
    for (k in seq_len(n_angles - 1)) {
      coords[k] <- sin_prod * cos(ang[k])
      sin_prod <- sin_prod * sin(ang[k])
    }
    coords[ncomp - 1] <- sin_prod * cos(ang[n_angles])
    coords[ncomp]     <- sin_prod * sin(ang[n_angles])
    child_corrected[i, ] <- r * coords
  }

  child_corrected <- sweep(child_corrected, 2, centroid_parent, "+")
  offset <- colMeans(child_corrected) - centroid_parent
  child_corrected <- sweep(child_corrected, 2, offset, "-")

  colnames(child_corrected) <- colnames(scores_child)
  rownames(child_corrected) <- rownames(scores_child)

  diagnostics <- data.frame(
    r_child, r_parent, r_hat, r_ratio, pred_ratio,
    radial_damp, angle_damp, blend_weight = loc_blend
  )

  message("✅ PCA correction complete!")

  # ============================================================
  # Plot results
  # ============================================================

  if (plot) {
    if (length(plot_dims) > 2) plot_dims <- plot_dims[1:2]
    pc1 <- plot_dims[1]; pc2 <- plot_dims[2]

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(4,4,3,1))
    plot(pca_parent$x[, pc1], pca_parent$x[, pc2],
         pch = 16, col = "firebrick", cex = 0.8,
         xlab = paste0("PC", pc1), ylab = paste0("PC", pc2),
         main = "Additive Spline PCA Correction")
    points(pca_child$x[, pc1], pca_child$x[, pc2],
           pch = 16, col = "darkgreen", cex = 0.8)
    points(child_corrected[, pc1], child_corrected[, pc2],
           pch = 16, col = "navy", cex = 0.8)
    legend("topright", legend = c("Parent", "Child", "Child Corrected"),
           col = c("firebrick", "darkgreen", "navy"), pch = 16)
  }

  invisible(list(
    child_corrected = child_corrected,
    diagnostics = diagnostics,
    models = list(ratio = model_ratio, angles = models_angle)
  ))
}
