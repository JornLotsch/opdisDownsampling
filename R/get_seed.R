#' Fast C++ Implementation of Seed Recovery
#'
#' @param range Optional vector of seeds to search
#' @param fallback_seed Seed to return if none found
#' @param max_search Maximum seed value to search
#' @param step_size Step size for range searching
#' @param batch_size Batch size for memory management
#' @param verbose Whether to print progress
#' @return Integer seed value
#' @export
get_seed_cpp <- function(range = NULL, fallback_seed = 42, max_search = 2147483647,
                         step_size = 50000, batch_size = 10000, verbose = TRUE) {

  if (!exists(".Random.seed", envir = globalenv())) {
    stop("No RNG state found.")
  }

  current <- get(".Random.seed", envir = globalenv())

  # Helper function using C++ backend
  searchRangeCpp <- function(search_values) {
    if (length(search_values) == 0) return(NULL)

    result <- findMatchingSeedsBatchCpp(
      as.integer(search_values),
      current,
      batchSize = batch_size,
      verbose = verbose
    )

    if (length(result) > 0) {
      return(result[1])
    }
    return(NULL)
  }

  # Search logic implementation
  if (is.null(range)) {
    if (verbose) {
      cat("Starting C++ seed search up to", format(max_search, scientific = FALSE), "...\n")
    }

    # Stage 1: Common seeds
    common_seeds <- c(1, 42, 123, 1234, 12345, 2023, 2024, 2025,
                      1:100, 1000, 5000, 9999, 99999)
    if (verbose) {
      cat("Checking common seeds (", length(common_seeds), " values)...\n")
    }

    result <- searchRangeCpp(common_seeds)
    if (!is.null(result)) {
      if (verbose) cat("Found seed:", result, "\n")
      return(result)
    }

    # Stage 2: Systematic search
    range_start <- 101
    while (range_start <= max_search) {
      range_end <- min(range_start + step_size - 1, max_search)

      if (verbose) {
        cat("Searching range:", format(range_start, scientific = FALSE),
            "to", format(range_end, scientific = FALSE),
            "(", format(range_end - range_start + 1, scientific = FALSE), " values)...\n")
      }

      result <- searchRangeCpp(range_start:range_end)
      if (!is.null(result)) {
        if (verbose) cat("Found seed:", result, "\n")
        return(result)
      }

      range_start <- range_end + 1
    }

    if (verbose) {
      cat("C++ search completed. No matching seed found within range.\n")
    }
  } else {
    # Use provided range
    if (verbose) {
      cat("Searching provided range (", length(range), " values)...\n")
    }
    result <- searchRangeCpp(range)
    if (!is.null(result)) {
      if (verbose) cat("Found seed:", result, "\n")
      return(result)
    }
  }

  if (verbose) {
    cat("Using fallback seed:", fallback_seed, "\n")
  }
  return(fallback_seed)
}

#' Enhanced get_seed function with C++ optimization
#' @param range Optional vector of seeds to search
#' @param fallback_seed Seed to return if none found
#' @param max_search Maximum seed value to search
#' @param step_size Step size for range searching
#' @param use_cpp Whether to use C++ implementation (default: TRUE)
#' @param ... Additional arguments passed to get_seed_cpp
#' @return Integer seed value
#' @export
get_seed <- function(range = NULL, fallback_seed = 42, max_search = 2147483647,
                     step_size = 50000, use_cpp = TRUE, ...) {

  if (use_cpp) {
    return(get_seed_cpp(range = range,
                        fallback_seed = fallback_seed,
                        max_search = max_search,
                        step_size = step_size,
                        ...))
  } else {
    # Fallback to original R implementation
    warning("Using slower R implementation")

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

      # Clean up large intermediate object immediately
      rm(seed_states)

      if (length(matches) > 0) return(search_values[matches[1]])
      return(NULL)
    }

    if (is.null(range)) {
      cat("Starting exhaustive seed search up to", format(max_search, scientific = FALSE), "...\n")

      # Stage 1: Common seeds (fast check)
      common_seeds <- c(1, 42, 123, 1234, 12345, 2023, 2024, 2025,
                        1:100, 1000, 5000, 9999, 99999)
      cat("Checking common seeds (", length(common_seeds), " values)...\n")
      result <- search_range(common_seeds)
      if (!is.null(result)) {
        cat("Found seed:", result, "\n")
        return(result)
      }

      # Stage 2: Dynamic expanding ranges
      range_start <- 101

      while (range_start <= max_search) {
        range_end <- min(range_start + step_size - 1, max_search)

        cat("Searching range:", format(range_start, scientific = FALSE),
            "to", format(range_end, scientific = FALSE),
            "(", format(range_end - range_start + 1, scientific = FALSE), " values)...\n")

        result <- search_range(range_start:range_end)
        if (!is.null(result)) {
          cat("Found seed:", result, "\n")
          return(result)
        }

        # Move to next range
        range_start <- range_end + 1
      }

      cat("Exhaustive search completed. No matching seed found within range 1 to",
          format(max_search, scientific = FALSE), "\n")

    } else {
      # Use provided range
      cat("Searching provided range (", length(range), " values)...\n")
      result <- search_range(range)
      if (!is.null(result)) {
        cat("Found seed:", result, "\n")
        return(result)
      }
    }

    cat("Using fallback seed:", fallback_seed, "\n")
    return(fallback_seed)
  }
}