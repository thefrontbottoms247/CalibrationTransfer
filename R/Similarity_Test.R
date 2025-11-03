#' @title Similarity_Test
#' @description
#' Compute multiple similarity and distance metrics between corresponding
#' columns (e.g., components or loadings) of two matrices A and B.
#' Optionally produce a barplot visualization for each metric across components.
#'
#' @param A Numeric matrix (e.g., parent loadings, reference spectra).
#' @param B Numeric matrix (e.g., child loadings, corrected spectra).
#' @param ncomp Number of components (columns) to compare.
#' @param plot Logical; if TRUE, displays a barplot of the similarity metrics.
#' @return
#' Data frame of similarity/distance metrics for each component,
#' with global RV coefficient stored as attribute.
#' @export
#' @examples
#' sim <- Similarity_Test(A = parent.loadings, B = child.loadings, ncomp = 5, plot = TRUE)

Similarity_Test <- function(A, B, ncomp = 5, plot = FALSE) {
  if (!is.matrix(A) || !is.matrix(B))
    stop("A and B must be matrices.")
  if (ncol(A) < ncomp || ncol(B) < ncomp)
    stop("Both matrices must have at least ncomp columns.")
  
  # --- Ensure numeric and aligned ---
  A <- as.matrix(A)
  B <- as.matrix(B)
  A <- A[, 1:ncomp, drop = FALSE]
  B <- B[, 1:ncomp, drop = FALSE]
  
  # --- Helper functions ---
  cosine <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
  pearson <- function(x, y) cor(x, y, method = "pearson")
  spearman <- function(x, y) cor(x, y, method = "spearman")
  euclidean <- function(x, y) sqrt(sum((x - y)^2))
  manhattan <- function(x, y) sum(abs(x - y))
  chebyshev <- function(x, y) max(abs(x - y))
  sam <- function(x, y) acos(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
  congruence <- function(x, y) sum(x * y) / sqrt(sum(x^2) * sum(y^2))
  
  # --- Compute metrics per component ---
  out <- data.frame(
    Component = 1:ncomp,
    Cosine = NA_real_,
    Pearson = NA_real_,
    Spearman = NA_real_,
    Euclidean = NA_real_,
    Manhattan = NA_real_,
    Chebyshev = NA_real_,
    SAM = NA_real_,
    Congruence = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:ncomp) {
    x <- A[, i]
    y <- B[, i]
    
    out[i, "Cosine"]     <- cosine(x, y)
    out[i, "Pearson"]    <- pearson(x, y)
    out[i, "Spearman"]   <- spearman(x, y)
    out[i, "Euclidean"]  <- euclidean(x, y)
    out[i, "Manhattan"]  <- manhattan(x, y)
    out[i, "Chebyshev"]  <- chebyshev(x, y)
    out[i, "SAM"]        <- sam(x, y)
    out[i, "Congruence"] <- congruence(x, y)
  }
  
  # --- Global RV coefficient (matrix-level correlation) ---
  if (requireNamespace("FactoMineR", quietly = TRUE)) {
    rv <- tryCatch(
      FactoMineR::RV(A[, 1:ncomp, drop = FALSE],
                     B[, 1:ncomp, drop = FALSE]),
      error = function(e) NA_real_
    )
  } else {
    rv <- NA_real_
  }
  
  attr(out, "RV") <- rv
  attr(out, "type") <- "Similarity_Test"
  
  # --- Plot option ---
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not installed — skipping plot.")
    } else {
      library(ggplot2)
      df_long <- reshape2::melt(out, id.vars = "Component")
      
      # Mark similarity vs distance metrics for scaling
      sim_metrics <- c("Cosine", "Pearson", "Spearman", "Congruence")
      df_long$Type <- ifelse(df_long$variable %in% sim_metrics,
                             "Similarity", "Distance")
      
      # Normalize distances so they can appear on same plot
      df_long$value[df_long$Type == "Distance"] <-
        scale(df_long$value[df_long$Type == "Distance"])
      
      ggplot(df_long, aes(x = factor(Component), y = value,
                          fill = variable)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
        labs(title = "Similarity and Distance Metrics per Component",
             subtitle = paste0("Global RV = ",
                               format(rv, digits = 3, nsmall = 3)),
             x = "Component", y = "Value") +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
        )
    }
  }
  
  return(out)
}
