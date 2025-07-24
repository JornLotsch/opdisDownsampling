#' @title Determine optimal number of cores for parallel processing
#' @description Calculates the number of CPU cores to use for parallel operations,
#' leaving one core free for system operations and respecting user-specified limits.
#' @param max_cores Integer or NULL. Maximum number of cores to use. If NULL,
#' uses all available cores minus one (default: NULL).
#' @return Integer. Number of cores to use for parallel processing (minimum 1).
#' @details This function ensures at least one core remains available for system
#' operations by using \code{parallel::detectCores() - 1}. If \code{max_cores} is
#' specified, it returns the minimum of the detected cores and the user limit.
#' @keywords internal
#' @importFrom parallel detectCores
#' @examples
#' \dontrun{
#' # Use all available cores minus one
#' n_cores <- determine_n_cores()
#'
#' # Limit to maximum 4 cores
#' n_cores <- determine_n_cores(max_cores = 4)
#' }
determine_n_cores <- function(max_cores = NULL) {
  if (!is.null(max_cores)) {
    return(max(1, min(parallel::detectCores() - 1, max_cores)))
  }
  max(1, parallel::detectCores() - 1)
}

#' @title Apply function with text progress bar
#' @description Sequential version of lapply that displays a text progress bar
#' to monitor computation progress, useful for long-running operations.
#' @param X List or vector. Input data to iterate over.
#' @param FUN Function. Function to apply to each element of X.
#' @param ... Additional arguments passed to FUN.
#' @return List. Results of applying FUN to each element of X.
#' @details Creates a text progress bar using \code{txtProgressBar} and updates
#' it after processing each element. The progress bar is automatically closed
#' upon completion. This function is particularly useful when parallel processing
#' is not available or desired.
#' @keywords internal
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' # Apply function with progress bar
#' results <- lapply_with_bar(1:100, function(x) {
#'   Sys.sleep(0.1)  # Simulate computation
#'   x^2
#' })
#' }
lapply_with_bar <- function(X, FUN, ...) {
  pb <- txtProgressBar(min = 0, max = length(X), style = 3)
  result <- vector("list", length(X))
  for (i in seq_along(X)) {
    result[[i]] <- FUN(X[[i]], ...)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  result
}

#' @title Find maximum value ignoring NaN/NA values
#' @description Finds the maximum value in data while ignoring NaN and NA values.
#' For matrices, returns column-wise maxima. For vectors, returns the overall maximum.
#' @param Data Numeric vector or matrix. Input data to find maximum values from.
#' @return Numeric. For vectors: single maximum value. For matrices: vector of
#' column-wise maximum values.
#' @details This function handles both vector and matrix inputs. For matrices,
#' it applies the maximum function column-wise using \code{apply} with
#' \code{na.rm = TRUE}. For vectors, it returns the single maximum value.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Vector example
#' vec <- c(1, 2, NA, 4, NaN, 3)
#' nanmax(vec)  # Returns 4
#'
#' # Matrix example
#' mat <- matrix(c(1, NA, 3, 2, NaN, 4), nrow = 2)
#' nanmax(mat)  # Returns column maxima
#' }
nanmax <- function(Data) {
  if (length(dim(Data)) == 2) {
    # Matrix case: return column-wise maxima
    SpaltenMinima <- apply(Data, 2, function(x) max(x, na.rm = TRUE))
    SpaltenInd <- NaN
  } else {
    # Vector case: return overall maximum
    SpaltenMinima <- max(Data, na.rm = TRUE)
    SpaltenInd <- which(Data == SpaltenMinima)
  }
  return(SpaltenMinima)
}

# ===== VALIDATION AND ERROR HANDLING FUNCTIONS =====

#' @title Validate reduced diagnostic results
#' @description Validates the structure and content of trial results from sample_and_analyze
#' @param ReducedDiag List of trial results to validate
#' @return Logical. TRUE if validation passes, otherwise stops with error
#' @keywords internal
validate_reduced_diag <- function(ReducedDiag) {
  if (length(ReducedDiag) == 0) {
    stop("opdisDownsampling: No valid results from sample_and_analyze.")
  }

  if (is.null(ReducedDiag[[1]])) {
    stop("opdisDownsampling: First trial result is NULL.")
  }

  if (!is.list(ReducedDiag[[1]]) ||
    is.null(ReducedDiag[[1]][[1]]) ||
    length(ReducedDiag[[1]][[1]]) == 0) {
    stop("opdisDownsampling: Invalid structure in trial results.")
  }

  TRUE
}

#' @title Validate statistical matrices
#' @description Validates statistical matrices for problematic values
#' @param matrices List of matrices to validate
#' @param matrix_names Character vector of matrix names for error messages
#' @return Invisible NULL. Issues warnings for problematic matrices
#' @keywords internal
validate_matrices <- function(matrices, matrix_names) {
  for (i in seq_along(matrices)) {
    if (all(is.na(matrices[[i]]))) {
      warning(sprintf("opdisDownsampling: Matrix '%s' contains all NA values.",
                      matrix_names[i]), call. = FALSE)
    }
  }
}

# ===== BEST TRIAL SELECTION FUNCTIONS =====

#' @title Select best trial for OptimizeBetween mode
#' @description Selects the best trial when optimizing between reduced and removed data
#' @param AD_reduced_vs_removed_statMat Matrix of statistics comparing reduced vs removed data
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_optimize_between <- function(AD_reduced_vs_removed_statMat) {
  if (is.null(AD_reduced_vs_removed_statMat) || nrow(AD_reduced_vs_removed_statMat) == 0) {
    stop("opdisDownsampling: No data available for OptimizeBetween selection.")
  }

  max_vals <- apply(AD_reduced_vs_removed_statMat, 1, max, na.rm = TRUE)
  if (all(is.infinite(max_vals))) {
    warning("opdisDownsampling: All values are infinite in reduced vs removed comparison.",
            call. = FALSE)
    return(1)
  }

  which.min(max_vals)
}

#' @title Select best trial for threefold optimization
#' @description Selects the best trial considering reduced, removed, and reduced-vs-removed comparisons
#' @param AD_reduced_statMat Matrix of statistics for reduced data vs original
#' @param AD_removed_statMat Matrix of statistics for removed data vs original
#' @param AD_reduced_vs_removed_statMat Matrix of statistics for reduced vs removed data
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_threefold <- function(AD_reduced_statMat, AD_removed_statMat,
                                        AD_reduced_vs_removed_statMat) {
  # Calculate ranks for all three matrices
  R_AD_reduced_statMat <- rank(apply(AD_reduced_statMat, 1, max, na.rm = TRUE),
                               ties.method = "first")
  R_AD_removed_statMat <- rank(apply(AD_removed_statMat, 1, max, na.rm = TRUE),
                               ties.method = "first")
  R_AD_reduced_vs_removed_statMat <- rank(apply(AD_reduced_vs_removed_statMat, 1, max, na.rm = TRUE),
                                          ties.method = "first")

  # Combine ranks
  R_AD_all_statMat <- cbind(R_AD_reduced_statMat, R_AD_removed_statMat, R_AD_reduced_vs_removed_statMat)

  # Find best trial
  best_trial <- which.min(apply(R_AD_all_statMat, 1, max, na.rm = TRUE))

  # Clean up temporary variables immediately
  rm(R_AD_reduced_statMat, R_AD_removed_statMat, R_AD_reduced_vs_removed_statMat, R_AD_all_statMat)

  return(best_trial)
}

#' @title Select best trial considering removed data
#' @description Selects the best trial considering both reduced and removed data vs original
#' @param AD_reduced_statMat Matrix of statistics for reduced data vs original
#' @param AD_removed_statMat Matrix of statistics for removed data vs original
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_check_removed <- function(AD_reduced_statMat, AD_removed_statMat) {
  # Calculate ranks for both matrices
  R_AD_reduced_statMat <- rank(apply(AD_reduced_statMat, 1, max, na.rm = TRUE),
                               ties.method = "first")
  R_AD_removed_statMat <- rank(apply(AD_removed_statMat, 1, max, na.rm = TRUE),
                               ties.method = "first")

  # Combine ranks
  R_AD_all_statMat <- cbind(R_AD_reduced_statMat, R_AD_removed_statMat)

  # Find best trial
  best_trial <- which.min(apply(R_AD_all_statMat, 1, max, na.rm = TRUE))

  # Clean up temporary variables immediately
  rm(R_AD_reduced_statMat, R_AD_removed_statMat, R_AD_all_statMat)

  return(best_trial)
}

#' @title Select best trial for reduced data only
#' @description Selects the best trial considering only reduced data vs original
#' @param AD_reduced_statMat Matrix of statistics for reduced data vs original
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_reduced_only <- function(AD_reduced_statMat) {
  if (is.null(AD_reduced_statMat) || nrow(AD_reduced_statMat) == 0) {
    stop("opdisDownsampling: No data available for reduced-only selection.")
  }

  max_vals <- apply(AD_reduced_statMat, 1, max, na.rm = TRUE)
  if (all(is.infinite(max_vals))) {
    warning("opdisDownsampling: All values are infinite in reduced matrix comparison.",
            call. = FALSE)
    return(1)
  }

  which.min(max_vals)
}