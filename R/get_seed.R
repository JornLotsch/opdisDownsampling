#' Get the Current Random Number Seed
#'
#' This function retrieves the current random number seed and the RNGkind for
#' reproducibility.
#'
#' @return A numeric vector containing the current random number seed.
#'
get_seed <- function() {
  # Check if the .Random.seed object exists in the global environment
  if (!exists(".Random.seed", envir = globalenv(), mode = "numeric", inherits = FALSE)) {
    # If not, generate a new random number and store it in .Random.seed
    runif(1L)
  }
  
  # Retrieve the current random number seed
  seed <- get(".Random.seed", envir = globalenv(), mode = "numeric", inherits = FALSE)
  
  # Attach the current RNGkind to the seed
  attr(seed, "RNGkind") <- RNGkind()
  
  return(seed)
}
