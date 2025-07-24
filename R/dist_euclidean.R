#' Calculate the Euclidean Distance of Two Probability Density Distributions
#'
#' This function calculates the Euclidean distance between two probability density
#' distributions expressed as Smooth Density Histograms.
#'
#' @param vector1 A numeric vector representing the first probability density distribution.
#' @param vector2 A numeric vector representing the second probability density distribution.
#'
#' @return The Euclidean distance between the two probability density distributions.
#'
#' @importFrom SmoothDensHist1dim SmoothDensHist1dim
#' @importFrom stats dist
#'
EucDist <- function(vector1, vector2) {
  # Create a Smooth Density Histogram for the combined input vectors
  densA <- SmoothDensHist1dim(c(vector1, vector2))
  xA <- densA$Kernels
  
  # Create Smooth Density Histograms for each input vector
  dens1 <- SmoothDensHist1dim(x = vector1, KernelsOrNbins = xA)
  dens2 <- SmoothDensHist1dim(x = vector2, KernelsOrNbins = xA)
  
  # Calculate the Euclidean distance between the two Smooth Density Histograms
  return(dist(rbind(dens1$SDH, dens2$SDH))[1])
}
