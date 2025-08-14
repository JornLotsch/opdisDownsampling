\name{get_seed}
\alias{get_seed}
\alias{get_seed_cpp}
\title{Random Seed Recovery from RNG State}
\description{
  Functions for recovering the original seed value that produced the current random number generator state. Provides both R and C++ implementations with the C++ version offering significantly improved performance for large search spaces.
}
\usage{
get_seed(range = NULL, fallback_seed = 42, max_search = 2147483647,
  step_size = 50000, use_cpp = TRUE, ...)

get_seed_cpp(range = NULL, fallback_seed = 42, max_search = 2147483647,
  step_size = 50000, batch_size = 10000, verbose = TRUE)
}
\arguments{
 \item{range}{Optional integer vector of specific seed values to search. If provided, only these seeds will be tested instead of systematic range searching.}

 \item{fallback_seed}{Integer seed value to return if no matching seed is found during the search process (default: 42).}

 \item{max_search}{Maximum seed value to search up to when performing systematic range searching. Must be a positive integer within the valid range for R's random number generator (default: 2147483647).}

 \item{step_size}{Step size for systematic range searching when no specific range is provided. Larger values speed up search but may miss the target seed if it falls between steps (default: 50000).}

 \item{use_cpp}{Logical; if \code{TRUE}, uses the fast C++ implementation via \code{get_seed_cpp()}. If \code{FALSE}, falls back to the slower R implementation with a warning (default: \code{TRUE}).}

 \item{batch_size}{Integer specifying the number of seeds to process in each C++ batch operation. Larger batches are more memory efficient but require more RAM. Only used in \code{get_seed_cpp()} (default: 10000).}

 \item{verbose}{Logical; if \code{TRUE}, prints progress information during the search process including batch progress and timing information. Only used in \code{get_seed_cpp()} (default: \code{TRUE}).}

 \item{...}{Additional arguments passed to \code{get_seed_cpp()} when \code{use_cpp = TRUE}.}
}

\value{
  Returns an integer representing the seed value that reproduces the current random number generator state.
  If no matching seed is found within the search parameters, returns the \code{fallback_seed} value.
}
\details{
  The functions work by systematically testing seed values to find one that reproduces the current RNG state stored in \code{.Random.seed}. The search process:

  \itemize{
    \item{Tests each candidate seed by setting it and comparing the resulting RNG state}
    \item{Uses efficient C++ implementation for faster processing of large search spaces}
    \item{Supports both targeted searching (via \code{range} parameter) and systematic range searching}
    \item{Employs batched processing to optimize memory usage and performance}
  }

  \strong{Performance Considerations:}

  The C++ implementation (\code{get_seed_cpp()}) provides significant performance improvements:
  \itemize{
    \item{Batch processing reduces overhead for large search spaces}
    \item{Optimized memory management prevents excessive RAM usage}
    \item{Native C++ random number generation matching R's implementation}
    \item{Progress reporting for long-running searches}
  }

  \strong{Search Strategy:}

  \itemize{
    \item{If \code{range} is provided: Tests only the specified seed values}
    \item{If \code{range} is NULL: Performs systematic search from 1 to \code{max_search} in steps of \code{step_size}}
    \item{Search terminates immediately when a matching seed is found}
    \item{Returns \code{fallback_seed} if no match is found within the search parameters}
  }

  \strong{Memory Management:}

  The C++ implementation uses batched processing controlled by \code{batch_size} to:
  \itemize{
    \item{Process large search ranges without excessive memory allocation}
    \item{Provide regular progress updates during long searches}
    \item{Allow interruption of long-running operations}
  }
}
\note{
  \itemize{
    \item{Requires an active RNG state (i.e., \code{.Random.seed} must exist)}
    \item{Large search ranges may take considerable time even with C++ optimization}
    \item{The search is deterministic but computationally intensive}
    \item{Consider using smaller \code{step_size} values if the initial search fails}
  }
}
\examples{
## Basic seed recovery after generating random numbers
set.seed(123)
recovered_seed <- get_seed()
print(recovered_seed)

}
\seealso{
  \code{\link{set.seed}}, \code{\link{.Random.seed}}
}
\keyword{random}
\keyword{seed}
\keyword{reproducibility}
\keyword{RNG}
\keyword{state recovery}