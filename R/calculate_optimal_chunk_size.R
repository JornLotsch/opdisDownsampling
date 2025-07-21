#' Calculate Optimal Chunk Size for Memory-Efficient Processing
#'
#' This function determines an optimal chunk size for processing trials based on
#' data characteristics, available memory, and computational resources.
#'
#' @param n_rows Number of rows in the dataset
#' @param n_cols Number of columns in the dataset (excluding class column)
#' @param nTrials Total number of trials to process
#' @param nProc Number of processor cores available
#'
#' @return An integer representing the optimal chunk size
#'
#' @details The function considers:
#' \itemize{
#'   \item Data size and estimated memory usage per trial
#'   \item Available system memory (on Linux systems)
#'   \item Number of processor cores
#'   \item Computational efficiency vs memory usage trade-offs
#' }
#'
#' @examples
#' \dontrun{
#' chunk_size <- calculate_optimal_chunk_size(1000, 50, 500, 4)
#' }
#'
#' @export
calculate_optimal_chunk_size <- function(n_rows, n_cols, nTrials, nProc) {
  # Calculate data characteristics
  data_size_mb <- (n_rows * n_cols * 8) / (1024^2)  # Approximate size in MB (8 bytes per double)

  # Estimate memory usage per trial
  # Each trial creates copies of subsets, so memory scales with data size
  memory_per_trial_mb <- data_size_mb * 1.5  # Factor for temporary objects and operations

  # Get available memory (rough estimate)
  # Try to use at most 25% of available memory for chunking
  available_mem_mb <- tryCatch({
    if (Sys.info()[["sysname"]] == "Linux") {
      # On Linux, try to read from /proc/meminfo
      meminfo <- readLines("/proc/meminfo", n = 3)
      available_line <- meminfo[grep("MemAvailable|MemFree", meminfo)[1]]
      if (length(available_line) > 0) {
        available_mem_kb <- as.numeric(sub(".*?([0-9]+).*", "\\1", available_line))
        available_mem_kb / 1024 * 0.25  # Use 25% of available memory
      } else {
        # Fallback if parsing fails
        max(1000, data_size_mb * 10)
      }
    } else {
      # Fallback for non-Linux systems: assume reasonable amount based on data size
      max(1000, data_size_mb * 10)  # At least 1GB or 10x data size
    }
  }, error = function(e) {
    # If memory detection fails, use conservative estimate
    max(1000, data_size_mb * 5)
  })

  # Calculate optimal chunk size based on memory constraints
  if (memory_per_trial_mb > 0) {
    max_parallel_trials <- max(1, floor(available_mem_mb / (memory_per_trial_mb * nProc)))
    chunk_size_memory <- min(nTrials, max_parallel_trials)
  } else {
    chunk_size_memory <- 50  # Default fallback
  }

  # Balance memory constraints with computational efficiency
  chunk_size <- if (nTrials <= 10) {
    nTrials  # Process all small jobs at once
  } else if (nTrials <= 50) {
    max(10, min(chunk_size_memory, nTrials))
  } else if (nTrials <= 500) {
    # Medium datasets: balance between memory and overhead
    max(20, min(chunk_size_memory, ceiling(nTrials / max(2, nProc))))
  } else {
    # Large datasets: prioritize memory efficiency
    if (data_size_mb > 500) {  # Large data (>500MB)
      max(10, min(25, chunk_size_memory))
    } else if (data_size_mb > 100) {  # Medium data (100-500MB)
      max(20, min(50, chunk_size_memory))
    } else {  # Smaller data (<100MB)
      max(50, min(100, chunk_size_memory))
    }
  }

  # Ensure chunk_size is reasonable
  chunk_size <- max(1, min(chunk_size, nTrials))

  return(chunk_size)
}

#' Print Chunk Size Diagnostics
#'
#' Helper function to print diagnostic information about chunk size calculation.
#' Useful for debugging and understanding memory usage patterns.
#'
#' @param n_rows Number of rows in the dataset
#' @param n_cols Number of columns in the dataset
#' @param nTrials Total number of trials
#' @param chunk_size Calculated chunk size
#' @param verbose Logical, whether to print diagnostics
#'
#' @return Invisible NULL
#'
#' @export
print_chunk_diagnostics <- function(n_rows, n_cols, nTrials, chunk_size, verbose = FALSE) {
  if (verbose) {
    data_size_mb <- (n_rows * n_cols * 8) / (1024^2)
    memory_per_trial_mb <- data_size_mb * 1.5
    n_chunks <- ceiling(nTrials / chunk_size)

    message(sprintf("Chunk size diagnostics:"))
    message(sprintf("  Data: %d rows x %d cols (%.1f MB)", n_rows, n_cols, data_size_mb))
    message(sprintf("  Estimated memory per trial: %.1f MB", memory_per_trial_mb))
    message(sprintf("  Trials: %d, Chunk size: %d, Number of chunks: %d",
                    nTrials, chunk_size, n_chunks))
  }
  invisible(NULL)
}