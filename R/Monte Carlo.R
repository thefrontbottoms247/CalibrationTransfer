#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export
Monte_Carlo <- function(parent, 
                       class, 
                       child_cor, 
                       R = 100, 
                       prop_train = 0.5, 
                       seed = 123) 
{
  set.seed(seed)
  library(pls)
  library(Metrics)
  library(ggplot2)
  
  # Convert to numeric and standardize column names
  parent <- as.data.frame(lapply(parent, as.numeric))
  child_cor <- as.data.frame(lapply(child_cor, as.numeric))
  class <- as.numeric(class)
  
  colnames_clean <- paste0("X", make.names(colnames(parent)))
  colnames(parent) <- colnames_clean
  colnames(child_cor) <- colnames_clean
  
  n <- nrow(parent)
  rmse_parent_self <- numeric(R)
  rmse_child_corr <- numeric(R)
  
  for (i in seq_len(R)) {
    idx_train <- sample(seq_len(n), size = floor(prop_train * n))
    idx_test <- setdiff(seq_len(n), idx_train)
    
    ncomp_use <- min(10, length(idx_train) - 1, ncol(parent))
    
    # Train model on parent
    train_data <- data.frame(y = class[idx_train], parent[idx_train, ])
    model <- plsr(y ~ ., data = train_data, ncomp = ncomp_use, validation = "none")
    
    # Test on parent (self performance)
    pred_parent <- drop(predict(model, newdata = parent[idx_test, ], ncomp = ncomp_use))
    rmse_parent_self[i] <- rmse(class[idx_test], pred_parent)
    
    # Test on child_cor (transfer performance)
    pred_child <- drop(predict(model, newdata = child_cor[idx_test, ], ncomp = ncomp_use))
    rmse_child_corr[i] <- rmse(class[idx_test], pred_child)
  }
  
  df.monte <<- data.frame(
    iter = seq_len(R),
    parent_self_rmse = rmse_parent_self,
    child_corr_rmse = rmse_child_corr,
    ratio = rmse_child_corr / rmse_parent_self
  )
  
  # Summary statistics
  cat("\n=== Parent Self-RMSE Summary ===\n")
  cat("Mean:   ", round(mean(rmse_parent_self), 4), "\n")
  cat("Median: ", round(median(rmse_parent_self), 4), "\n")
  cat("SD:     ", round(sd(rmse_parent_self), 4), "\n")
  cat("Min:    ", round(min(rmse_parent_self), 4), "\n")
  cat("Max:    ", round(max(rmse_parent_self), 4), "\n")
  
  cat("\n=== Percentile Ranges ===\n")
  quantiles <- quantile(rmse_parent_self, probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1))
  for (i in seq_along(quantiles)) {
    cat(sprintf("%3.0f%%: %7.4f\n", as.numeric(sub("%", "", names(quantiles)[i])), quantiles[i]))
  }
  
  cat("\n=== Percentage within ranges ===\n")
  ranges <- list(c(0, 2), c(2, 2.5), c(2.5, 3), c(3, 3.5), c(3.5, 4), c(4, Inf))
  for (r in ranges) {
    count <- sum(rmse_parent_self >= r[1] & rmse_parent_self < r[2])
    pct <- count / R * 100
    cat(sprintf("[%.1f, %.1f): %3d iterations (%5.1f%%)\n", r[1], r[2], count, pct))
  }
  
  # Violin plot
  plot_data <- data.frame(
    RMSE = c(rmse_parent_self, rmse_child_corr),
    Type = factor(rep(c("Parent Self", "Child Corrected"), each = R),
                  levels = c("Parent Self", "Child Corrected"))
  )
  
  p <- ggplot(plot_data, aes(x = Type, y = RMSE, fill = Type)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
    labs(title = "Monte Carlo RMSE Distribution",
         subtitle = paste0("R = ", R, " iterations"),
         x = NULL, y = "RMSE") +
    theme_minimal() +
    theme(legend.position = "none", text = element_text(size = 12)) +
    scale_fill_manual(values = c("Parent Self" = "#3498db", "Child Corrected" = "#e74c3c"))
  
  print(p)
}