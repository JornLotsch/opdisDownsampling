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

#' @title Select best trial for multiple matrices
#' @description Selects the best trial considering any number of comparison matrices
#' @param ... Matrices of statistics for various comparisons (e.g., reduced vs original, removed vs original, etc.)
#' @param matrix_names Optional character vector of matrix names for error messages.
#'   If not provided, will use argument names or generate default names.
#' @param WorstSample A logical value for testing purpose reversing the split ranking to obtain
#'   the least similar subsample (default: FALSE).
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_multi <- function(..., matrix_names = NULL, WorstSample = FALSE) {

  # Collect all matrices
  matrices <- list(...)

  # Validate we have at least one matrix
  if (length(matrices) == 0) {
    stop("opdisDownsampling: No matrices provided for selection.")
  }

  # Generate matrix names if not provided
  if (is.null(matrix_names)) {
    arg_names <- match.call(expand.dots = FALSE)$`...`
    if (!is.null(arg_names)) {
      matrix_names <- sapply(arg_names, deparse)
    } else {
      matrix_names <- paste0("Matrix_", seq_along(matrices))
    }
  }

  # Validate matrix names length matches matrices
  if (length(matrix_names) != length(matrices)) {
    matrix_names <- paste0("Matrix_", seq_along(matrices))
  }

  # Input validity checks
  for (i in seq_along(matrices)) {
    if (is.null(matrices[[i]]) || nrow(matrices[[i]]) == 0) {
      stop(paste0("opdisDownsampling: No data available for multi-matrix selection in ",
                  matrix_names[i], "."))
    }
  }

  # Calculate ranks for all matrices
  rank_results <- lapply(matrices, function(mat) {
    rank(apply(mat, 1, max, na.rm = TRUE), ties.method = "first")
  })

  # Combine ranks into a matrix
  R_all_statMat <- do.call(cbind, rank_results)

  # Find best trial
  best_trial <- if (!WorstSample) {
    which.min(apply(R_all_statMat, 1, max, na.rm = TRUE))
  } else {
    which.max(apply(R_all_statMat, 1, max, na.rm = TRUE))
  }

  return(best_trial)
}