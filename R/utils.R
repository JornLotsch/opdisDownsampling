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