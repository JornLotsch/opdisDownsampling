#' Calculate the Euclidean Distance of Two Probability Density Distributions
#'
#' This internal helper calculates the Euclidean distance between two probability
#' density distributions represented as smooth density histograms.
#'
#' @param vector1 A numeric vector representing the first distribution.
#' @param vector2 A numeric vector representing the second distribution.
#'
#' @return The Euclidean distance between the two smoothed density histograms.
#'
#' @importFrom stats dist
#'
EucDist <- function(vector1, vector2) {
  # Create a Smooth Density Histogram for the combined input vectors
  densA <- SmoothDensHist1dim(c(vector1, vector2))
  xA <- densA$Kernels

  # Create Smooth Density Histograms for each input vector
  dens1 <- SmoothDensHist1dim(x = vector1, KernelsOrNbins = xA)
  dens2 <- SmoothDensHist1dim(x = vector2, KernelsOrNbins = xA)

  # Merge the Smooth Density Histograms into a single data frame
  dfPDF <- data.frame(
    Kernels = xA,
    SDH.x = dens1$SDH,
    SDH.y = dens2$SDH
  )

  # Calculate and return the Euclidean distance between the two smoothed densities
  return(as.numeric(stats::dist(dfPDF[, c("SDH.x", "SDH.y")])))
}
