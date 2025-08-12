#' Get the Current Random Number Seed by Reverse Engineering
#'
#' This function attempts to recover the seed value that was used with \code{set.seed()}
#' to generate the current random number generator (RNG) state. It works by testing
#' seed values to find which one produces the same .Random.seed state as the current one.
#' The function uses a smart search strategy that tests common seed values first before
#' performing a systematic search.
#'
#' @param range A numeric vector specifying the range of seed values to test.
#'   Default is NULL, which uses a smart search strategy: tests common seed values
#'   (1, 42, 123, etc.) first, then performs systematic search from 101 to 10000.
#'   You can provide a custom range for specific search requirements.
#' @param fallback_seed An integer value to return if no matching seed is found
#'   within the specified range. Default is 42.
#' @param verbose A logical value indicating whether to print search progress
#'   and diagnostic information. Default is FALSE.
#'
#' @return An integer representing the seed value that produces the current RNG state,
#'   or the fallback value if no matching seed is found within the specified range.
#'
#' @details
#' The function works by:
#' \itemize{
#'   \item Checking if .Random.seed exists in the global environment
#'   \item Storing the current .Random.seed state
#'   \item Using a smart search strategy:
#'     \itemize{
#'       \item First testing common seed values (1, 42, 123, 1234, etc.)
#'       \item Then performing systematic search in a reasonable range (101-10000)
#'     }
#'   \item For each value, calling \code{set.seed()} with that value
#'   \item Comparing the resulting .Random.seed with the saved current state
#'   \item Returning the first matching seed value found, or the fallback value
#' }
#'
#' \strong{Performance Notes:}
#' \itemize{
#'   \item The smart search strategy significantly improves performance for common use cases
#'   \item Most commonly used seeds (1, 42, 123, etc.) are found very quickly
#'   \item The default search range is much smaller than the original (10,000 vs 100,000,000)
#'   \item For larger search ranges, consider using verbose = TRUE to monitor progress
#' }
#'
#' \strong{Important Notes:}
#' \itemize{
#'   \item This function requires that .Random.seed exists in the global environment
#'   \item The search range should be chosen based on the likely seed values used
#'   \item The function modifies .Random.seed during the search process
#'   \item This is a brute-force reverse engineering approach
#'   \item Always returns a valid integer, never NA
#' }
#'
#' @seealso \code{\link{set.seed}}
#'
#' @export
get_seed <- function(range = NULL, fallback_seed = 42, verbose = FALSE) {
  if (!exists(".Random.seed", envir = globalenv())) {
    stop("No RNG state found.")
  }

  current <- get(".Random.seed", envir = globalenv())

  # Define smart search order: common values first, then systematic search
  if (is.null(range)) {
    # Common seed values people use
    common_seeds <- c(1, 42, 123, 1234, 12345, 2023, 2024, 2025,
                      1:100, 1000, 5000, 9999, 99999)
    # Then systematic search in smaller range
    systematic_range <- 101:10000
    range <- c(common_seeds, systematic_range)
    # Remove duplicates while preserving order
    range <- range[!duplicated(range)]
  }

  total_tests <- length(range)
  if (verbose) {
    cat(sprintf("Searching %d seed values...\n", total_tests))
  }

  # Search with progress updates
  for (i in seq_along(range)) {
    s <- range[i]
    set.seed(s)

    if (identical(.Random.seed, current)) {
      if (verbose) {
        cat(sprintf("Found seed %d after %d tests\n", s, i))
      }
      return(s)
    }

    # Progress indication for long searches
    if (verbose && i %% 5000 == 0) {
      cat(sprintf("Tested %d/%d seeds...\n", i, total_tests))
    }
  }

  if (verbose) {
    cat(sprintf("No matching seed found in %d tests, returning fallback: %d\n",
                total_tests, fallback_seed))
  }

  return(fallback_seed)
}