#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export
MSID <- function(X, k, metric = "euclidean", max_iter = 50, seed = 123) {
  set.seed(seed)
  Xs <- scale(X)  # standardize for fair distance computation
  D <- as.matrix(dist(Xs, method = metric))
  n <- nrow(D)
  
  best_subset <- NULL
  best_min_dist <- -Inf
  
  for (iter in seq_len(max_iter)) {
    # random start
    subset <- sample(seq_len(n), k)
    improved <- TRUE
    iter_count <- 0
    
    while (improved && iter_count < 1000) {
      improved <- FALSE
      iter_count <- iter_count + 1
      
      # compute smallest pairwise distance among current subset
      current_min <- min(D[subset, subset][lower.tri(D[subset, subset])])
      
      # try swapping each subset element with a random outside element
      for (i in subset) {
        outsiders <- setdiff(seq_len(n), subset)
        for (o in outsiders) {
          candidate <- c(setdiff(subset, i), o)
          cand_min <- min(D[candidate, candidate][lower.tri(D[candidate, candidate])])
          if (cand_min > current_min) {
            subset <- candidate
            current_min <- cand_min
            improved <- TRUE
            break
          }
        }
        if (improved) break
      }
    }
    
    # retain best subset found
    if (current_min > best_min_dist) {
      best_min_dist <- current_min
      best_subset <- subset
    }
  }
  
  cat(sprintf("Best min interpoint distance = %.3f\n", best_min_dist))
  return(sort(best_subset))
}
