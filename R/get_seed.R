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
#'   (1, 42, 123, etc.) first, then performs systematic search from 101 to 100000.
#'   You can provide a custom range for specific search requirements.
#' @param fallback_seed An integer value to return if no matching seed is found
#'   within the specified range. Default is 42.
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
#'       \item Then performing systematic search from 101 to 100000
#'     }
#'   \item For each range, calling \code{set.seed()} with each value in that range
#'   \item Comparing the resulting .Random.seed with the saved current state
#'   \item Returning the first matching seed value found, or the fallback value
#' }
#'
#' \strong{Performance Notes:}
#' \itemize{
#'   \item The smart search strategy significantly improves performance for common use cases
#'   \item Most commonly used seeds (1, 42, 123, etc.) are found very quickly
#'   \item Uses efficient helper function to eliminate code duplication
#'   \item Always returns a valid integer, never NA
#' }
#'
#' \strong{Important Notes:}
#' \itemize{
#'   \item This function requires that .Random.seed exists in the global environment
#'   \item The search range should be chosen based on the likely seed values used
#'   \item The function modifies .Random.seed during the search process
#'   \item This is a brute-force reverse engineering approach
#'   \item Memory and time intensive for large search ranges
#' }
#'
#' @seealso \code{\link{set.seed}}
#'
#' @export
get_seed <- function(range = NULL, fallback_seed = 42) {
  if (!exists(".Random.seed", envir = globalenv())) {
    stop("No RNG state found.")
  }

  current <- get(".Random.seed", envir = globalenv())

  # Helper function to search for matching seed in a given range
  search_range <- function(search_values) {
    seed_states <- lapply(search_values, function(s) {
      set.seed(s)
      .Random.seed
    })
    matches <- which(vapply(seed_states, function(x) identical(x, current), logical(1)))
    if (length(matches) > 0) return(search_values[matches[1]])
    return(NULL)
  }

  # Define search ranges
  if (is.null(range)) {
    # Smart search order: common values first, then systematic search
    common_seeds <- c(1, 42, 123, 1234, 12345, 2023, 2024, 2025,
                      1:100, 1000, 5000, 9999, 99999)
    systematic_range <- 101:100000

    # Search common seeds first
    result <- search_range(common_seeds)
    if (!is.null(result)) return(result)

    # Then search systematic range
    result <- search_range(systematic_range)
    if (!is.null(result)) return(result)
  } else {
    # Use provided range
    result <- search_range(range)
    if (!is.null(result)) return(result)
  }

  return(fallback_seed)
}
