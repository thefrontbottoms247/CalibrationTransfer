#' @title Monte Carlo PLS Transfer Evaluation (Extended)
#' @description
#' Performs repeated random sub-sampling validation comparing:
#' Parent Self, Child Uncorrected, and Child Corrected spectra using PLS models.
#' Reports RMSE distributions, percent differences, and visual summaries.
#'
#' @param parent Matrix or data frame of parent features.
#' @param class Numeric or factor vector of response values or class labels.
#' @param child_uncor Matrix or data frame of uncorrected child features.
#' @param child_cor Matrix or data frame of corrected child features.
#' @param R Integer. Number of Monte Carlo iterations (default = 100).
#' @param prop_train Proportion of samples to use for training (default = 0.5).
#' @param seed Random seed for reproducibility (default = 123).
#' @param ncomp Number of PLS components to use (default = 10).
#'
#' @return A data frame (`df.monte`) summarizing RMSE distributions and percent differences.
#' @export
Monte_Carlo <- function(parent,
                        class,
                        child_uncor,
                        child_cor,
                        R = 100,
                        prop_train = 0.5,
                        seed = 123,
                        ncomp = 10) {

  set.seed(seed)
  library(pls)
  library(Metrics)
  library(ggplot2)

  # --- Input prep ---
  to_df <- function(x) {
    if (is.matrix(x)) x <- as.data.frame(x)
    as.data.frame(lapply(x, as.numeric))
  }

  parent     <- to_df(parent)
  child_uncor <- to_df(child_uncor)
  child_cor  <- to_df(child_cor)
  class      <- as.numeric(class)

  if (nrow(parent) != length(class) ||
      nrow(child_uncor) != length(class) ||
      nrow(child_cor) != length(class))
    stop("Dimension mismatch between predictors and class vector.")

  colnames_clean <- paste0("X", make.names(colnames(parent)))
  colnames(parent) <- colnames_clean
  colnames(child_uncor) <- colnames_clean
  colnames(child_cor) <- colnames_clean

  n <- nrow(parent)
  rmse_parent_self <- numeric(R)
  rmse_child_uncor <- numeric(R)
  rmse_child_cor   <- numeric(R)
  ncomp_actual     <- integer(R)

  cat("Dataset loaded: n =", n, "samples, p =", ncol(parent), "features\n")

  for (i in seq_len(R)) {
    idx_train <- sample(seq_len(n), size = floor(prop_train * n))
    idx_test  <- setdiff(seq_len(n), idx_train)
    ncomp_use <- min(ncomp, ncol(parent), length(idx_train) - 1)
    ncomp_use <- max(1, ncomp_use)
    ncomp_actual[i] <- ncomp_use

    train_data <- data.frame(y = class[idx_train], parent[idx_train, ])
    model <- plsr(y ~ ., data = train_data, ncomp = ncomp_use, validation = "none")

    pred_parent <- drop(predict(model, newdata = parent[idx_test, ], ncomp = ncomp_use))
    pred_uncor  <- drop(predict(model, newdata = child_uncor[idx_test, ], ncomp = ncomp_use))
    pred_cor    <- drop(predict(model, newdata = child_cor[idx_test, ], ncomp = ncomp_use))

    rmse_parent_self[i] <- rmse(class[idx_test], pred_parent)
    rmse_child_uncor[i] <- rmse(class[idx_test], pred_uncor)
    rmse_child_cor[i]   <- rmse(class[idx_test], pred_cor)
  }

  # --- Summaries ---
  df.monte <<- data.frame(
    iter = seq_len(R),
    parent_self_rmse = rmse_parent_self,
    child_uncor_rmse = rmse_child_uncor,
    child_corr_rmse  = rmse_child_cor,
    ratio_uncor = rmse_child_uncor / rmse_parent_self,
    ratio_cor   = rmse_child_cor / rmse_parent_self
  )

  df.monte$delta_uncor_pct <<- 100 * (rmse_child_uncor - rmse_parent_self) / rmse_parent_self
  df.monte$delta_cor_pct   <<- 100 * (rmse_child_cor - rmse_parent_self) / rmse_parent_self

  cat("\n=== Parent Self-RMSE Summary ===\n")
  cat("Mean:", round(mean(rmse_parent_self), 4),
      "| SD:", round(sd(rmse_parent_self), 4),
      "| Median:", round(median(rmse_parent_self), 4),
      "| Range:", paste0(round(range(rmse_parent_self), 4), collapse = " - "), "\n")

  cat("\n=== Percent difference (Uncorrected Child vs Parent) ===\n")
  delta_u <- df.monte$delta_uncor_pct
  cat("Mean Δ%:", round(mean(delta_u), 3),
      "| Median Δ%:", round(median(delta_u), 3),
      "| SD Δ%:", round(sd(delta_u), 3),
      "| Range:", paste0(round(range(delta_u), 3), collapse = " - "), "\n")

  cat("\n=== Percent difference (Corrected Child vs Parent) ===\n")
  delta_c <- df.monte$delta_cor_pct
  cat("Mean Δ%:", round(mean(delta_c), 3),
      "| Median Δ%:", round(median(delta_c), 3),
      "| SD Δ%:", round(sd(delta_c), 3),
      "| Range:", paste0(round(range(delta_c), 3), collapse = " - "), "\n")

  thr <- c(0.25, 0.5, 1, 2, 5)
  cat("\n=== Share within ±x% (Corrected) ===\n")
  for (t in thr) {
    pct <- mean(abs(delta_c) <= t) * 100
    cat(sprintf("±%-5.2f%% : %5.1f%% of iterations\n", t, pct))
  }
  cat("\nImproved (Child Corrected < Parent):", sum(delta_c < 0),
      "| Worse (Child Corrected > Parent):", sum(delta_c > 0), "\n")

  # --- Violin plot ---
  plot_data <- data.frame(
    RMSE = c(rmse_parent_self, rmse_child_uncor, rmse_child_cor),
    Type = factor(rep(c("Parent Self", "Child Uncorrected", "Child Corrected"), each = R),
                  levels = c("Parent Self", "Child Uncorrected", "Child Corrected"))
  )

  actual_lab <- if (length(unique(ncomp_actual)) == 1)
    as.character(unique(ncomp_actual))
  else paste0(median(ncomp_actual), " (median actual)")

  p <- ggplot(plot_data, aes(x = Type, y = RMSE, fill = Type)) +
    geom_violin(trim = FALSE, alpha = 0.8, color = "black", width = 1.1) +
    geom_boxplot(width = 0.12, fill = "white", alpha = 0.7, outlier.size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, fill = "red") +
    labs(
      title = "Monte Carlo RMSE Distribution",
      subtitle = paste0("R = ", R, " iterations, ncomp = ", actual_lab),
      x = NULL, y = "RMSE"
    ) +
    scale_fill_manual(values = c("Parent Self" = "#4DA3FF",
                                 "Child Uncorrected" = "#FFD580",
                                 "Child Corrected" = "#F08080")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  print(p)

  invisible(df.monte)
}
