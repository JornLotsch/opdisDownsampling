# ===== BEST TRIAL SELECTION FUNCTIONS =====

#' @title Select best trial for one matrix
#' @description Selects the best trial when optimizing for one matrix
#'   such as for the reduced data subset versus the original or for 
#'   reduced versus removed data
#' @param AD_statMat Matrix of statistics comparing data
#' @param WorstSample A logical value for testing purpose reversing the split ranking to obtain
#'   the least similar subsample (default: FALSE).
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_one_matrix <- function(AD_statMat, WorstSample = FALSE) {

  if (is.null(AD_statMat) || nrow(AD_statMat) == 0) {
    stop("opdisDownsampling: No data available for OptimizeBetween selection.")
  }

  max_vals <- apply(AD_statMat, 1, if (!WorstSample) max else min, na.rm = TRUE)

  if (all(is.infinite(max_vals))) {
    warning("opdisDownsampling: All values are infinite in reduced vs removed comparison.",
            call. = FALSE)
    return(1)
  }

  if (!WorstSample) which.min(max_vals) else which.max(max_vals)
}

#' @title Select best trial for threefold optimization
#' @description Selects the best trial considering reduced, removed, and reduced-vs-removed comparisons
#' @param AD_reduced_statMat Matrix of statistics for reduced data vs original
#' @param AD_removed_statMat Matrix of statistics for removed data vs original
#' @param AD_reduced_vs_removed_statMat Matrix of statistics for reduced vs removed data
#' @param WorstSample A logical value for testing purpose reversing the split ranking to obtain
#'   the least similar subsample (default: FALSE).
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_threefold <- function(AD_reduced_statMat, AD_removed_statMat,
                                        AD_reduced_vs_removed_statMat, WorstSample = FALSE) {

  # Input validity checks
  matrices <- list(AD_reduced_statMat, AD_removed_statMat, AD_reduced_vs_removed_statMat)
  matrix_names <- c("AD_reduced_statMat", "AD_removed_statMat", "AD_reduced_vs_removed_statMat")

  for (i in seq_along(matrices)) {
    if (is.null(matrices[[i]]) || nrow(matrices[[i]]) == 0) {
      stop(paste0("opdisDownsampling: No data available for threefold selection in ", matrix_names[i], "."))
    }
  }

  # Calculate ranks for all three matrices
  rank_results <- lapply(matrices, function(mat) rank(apply(mat, 1, max, na.rm = TRUE), ties.method = "first"))
  R_AD_reduced_statMat <- rank_results[[1]]
  R_AD_removed_statMat <- rank_results[[2]]
  R_AD_reduced_vs_removed_statMat <- rank_results[[3]]

  # Combine ranks
  R_AD_all_statMat <- cbind(R_AD_reduced_statMat, R_AD_removed_statMat, R_AD_reduced_vs_removed_statMat)

  # Find best trial
  best_trial <- if (!WorstSample) which.min(apply(R_AD_all_statMat, 1, if (!WorstSample) max else min, na.rm = TRUE)) else 
    which.max(apply(R_AD_all_statMat, 1, max, na.rm = TRUE))

  # Clean up temporary variables immediately
  rm(R_AD_reduced_statMat, R_AD_removed_statMat, R_AD_reduced_vs_removed_statMat, R_AD_all_statMat)

  return(best_trial)
}

#' @title Select best trial considering removed data
#' @description Selects the best trial considering both reduced and removed data vs original
#' @param AD_reduced_statMat Matrix of statistics for reduced data vs original
#' @param AD_removed_statMat Matrix of statistics for removed data vs original
#' @param WorstSample A logical value for testing purpose reversing the split ranking to obtain
#'   the least similar subsample (default: FALSE).
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_check_removed <- function(AD_reduced_statMat, AD_removed_statMat, WorstSample = FALSE) {

  # Input validity checks
  matrices <- list(AD_reduced_statMat, AD_removed_statMat)
  matrix_names <- c("AD_reduced_statMat", "AD_removed_statMat")

  for (i in seq_along(matrices)) {
    if (is.null(matrices[[i]]) || nrow(matrices[[i]]) == 0) {
      stop(paste0("opdisDownsampling: No data available for check removed selection in ", matrix_names[i], "."))
    }
  }

  # Calculate ranks for both matrices
  rank_results <- lapply(matrices, function(mat) rank(apply(mat, 1, if (!WorstSample) max else min, na.rm = TRUE), ties.method = "first"))
  R_AD_reduced_statMat <- rank_results[[1]]
  R_AD_removed_statMat <- rank_results[[2]]

  # Combine ranks
  R_AD_all_statMat <- cbind(R_AD_reduced_statMat, R_AD_removed_statMat)

  # Find best trial
  best_trial <- if (!WorstSample) which.min(apply(R_AD_all_statMat, 1, max, na.rm = TRUE)) else 
    which.max(apply(R_AD_all_statMat, 1, max, na.rm = TRUE))
  
  # Clean up temporary variables immediately
  rm(R_AD_reduced_statMat, R_AD_removed_statMat, R_AD_all_statMat)

  return(best_trial)
}
