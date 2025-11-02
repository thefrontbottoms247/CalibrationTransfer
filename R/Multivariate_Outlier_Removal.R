#' @title Multivariate Outlier Removal
#' @description Detects and flags spectral outliers using multiple statistical and machine learning methods (e.g., Mahalanobis, Hotelling’s T², Q-residuals, k-means, DBSCAN, Isolation Forest, LOF, etc.).
#' @param spectra Numeric matrix or data frame containing spectral data (rows = samples, columns = variables).
#' @param method Character; specific method name ("mahalanobis", "euclidean", "hotelling", "qresiduals", "robust", "kmeans", "dbscan", "iso", "lof") or "all" to run all available methods.
#' @return Invisibly returns nothing; detected outlier indices are assigned to global variables (`outliers.method` for each method).
#' @export
Multivariate_Outlier_Removal <- function(spectra, method = "all") {
  
  library(solitude)   # Isolation Forest
  library(dbscan)     # for LOF
  library(keras)      # for autoencoder (optional, heavy)
  library(stats)
  library(mvtnorm)
  library(MASS)
  library(rrcov)
  library(cluster)
  library(dbscan)
  library(isotree)
  
  # Initialize progress bar if method is "all"
  if (method == "all") {
    pb <- txtProgressBar(min = 0, max = 9, style = 3)
    progress_counter <- 0
  }
  
  if (method == "mahalanobis" || method == "all") {
    # Mahalanobis Distance
    pca <- prcomp(spectra, scale. = TRUE)
    scores <- pca$x[, 1:10]  # top 10 PCs
    md <- mahalanobis(scores, colMeans(scores), cov(scores))
    cutoff <- qchisq(0.975, df = ncol(scores))
    outliers.mahalanobis <<- which(md > cutoff)
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "mahalanobis") return(invisible())
  }
  
  if (method == "euclidean" || method == "all") {
    # Euclidean Distance from centroid
    centroid <- colMeans(spectra)
    ed <- apply(spectra, 1, function(x) sqrt(sum((x - centroid)^2)))
    threshold <- mean(ed) + 3*sd(ed)
    outliers.euclidean <<- which(ed > threshold)
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "euclidean") return(invisible())
  }
  
  if (method == "hotelling" || method == "all") {
    # Hotelling's T²
    pca <- prcomp(spectra, scale. = TRUE)
    scores <- pca$x[, 1:10]
    cov_scores <- cov(scores[, 1:5])  # use first 5 PCs
    t2 <- mahalanobis(scores[, 1:5], colMeans(scores[, 1:5]), cov_scores)
    t2_cutoff <- qchisq(0.975, df = 5)
    outliers.hotelling <<- which(t2 > t2_cutoff)
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "hotelling") return(invisible())
  }
  
  if (method == "qresiduals" || method == "all") {
    # Q-residuals
    residuals <- sweep(spectra, 2, colMeans(spectra))
    qres <- rowSums(residuals^2)
    q_cutoff <- mean(qres) + 3*sd(qres)
    outliers.qresiduals <<- which(qres > q_cutoff)
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "qresiduals") return(invisible())
  }
  
  if (method == "robust" || method == "all") {
    # Robust Mahalanobis Distance using MCD
    pca <- prcomp(spectra, scale. = TRUE)
    rob <- CovMcd(pca$x[, 1:10])
    md_rob <- mahalanobis(pca$x[, 1:10], rob@center, rob@cov)
    cutoff_rob <- qchisq(0.975, df = ncol(pca$x[, 1:10]))
    outliers.robust <<- which(md_rob > cutoff_rob)
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "robust") return(invisible())
  }
  
  if (method == "kmeans" || method == "all") {
    # k-means distance
    kmeans_model <- kmeans(spectra, centers = 3, nstart = 10)
    dist_to_center <- sqrt(rowSums((spectra - kmeans_model$centers[kmeans_model$cluster, ])^2))
    thresh_kmeans <- mean(dist_to_center) + 3*sd(dist_to_center)
    outliers.kmeans <<- which(dist_to_center > thresh_kmeans)
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "kmeans") return(invisible())
  }
  
  if (method == "dbscan" || method == "all") {
    # DBSCAN
    db <- dbscan(spectra, eps = 10, minPts = 5) # eps may need tuning
    outliers.dbscan <<- which(db$cluster == 0)   # noise points
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "dbscan") return(invisible())
  }
  
  if (method == "iso" || method == "all") {
    # Isolation Forest
    pca <- prcomp(spectra, scale. = TRUE)
    pca_scores <- pca$x[, 1:10]
    
    # Set sample_size based on data size
    n_rows <- nrow(pca_scores)
    sample_size <- min(256, n_rows)
    
    iso <- isolationForest$new(sample_size = sample_size)
    iso$fit(pca_scores)
    pred_iso <- iso$predict(pca_scores)
    outliers.iso <<- which(pred_iso$anomaly_score > quantile(pred_iso$anomaly_score, 0.975))
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    if (method == "iso") return(invisible())
  }
  
  if (method == "lof" || method == "all") {
    # Local Outlier Factor (LOF)
    lof_scores <- lof(spectra, minPts = 10)
    thresh_lof <- mean(lof_scores) + 3*sd(lof_scores)
    outliers.lof <<- which(lof_scores > thresh_lof)
    if (method == "all") {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
      close(pb)
    }
    if (method == "lof") return(invisible())
  }
  
  return(invisible())
}