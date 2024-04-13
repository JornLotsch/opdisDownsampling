#' Calculate the Mean Relative Difference for Two Probability Density Distributions
#'
#' This function calculates the mean relative difference for two probability density
#' distributions expressed as Smooth Density Histograms.
#'
#' @param vector1 A numeric vector representing the first probability density distribution.
#' @param vector2 A numeric vector representing the second probability density distribution.
#'
#' @return The mean relative difference between the two probability density distributions.
#'
#' @importFrom SmoothDensHist1dim SmoothDensHist1dim
#'
amrdd <- function(vector1, vector2) {
  # Combine the two input vectors and create a Smooth Density Histogram
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
  
  # Calculate the absolute difference and the mean of the two distributions
  dfPDF$Diff <- abs(dfPDF$SDH.x - dfPDF$SDH.y)
  dfPDF$Mean <- (dfPDF$SDH.x + dfPDF$SDH.y) / 2
  dfPDF$Mean[dfPDF$Mean == 0] <- 1
  
  # Calculate the relative difference and return the mean
  dfPDF$RelDiff <- dfPDF$Diff / dfPDF$Mean
  return(mean(dfPDF$RelDiff, na.rm = TRUE))
}
