#Calculates the Euclidean distance of two  probability density distributions
#expressed as Smooth Density Histogram
#' @importFrom SmoothDensHist1dim SmoothDensHist1dim
EucDist <- function(vector1, vector2) {
  densA <- SmoothDensHist1dim(c(vector1, vector2))
  xA = densA$Kernels
  dens1 <- SmoothDensHist1dim(x = vector1, KernelsOrNbins = xA)
  y1 = data.frame(dens1)
  dens2 <- SmoothDensHist1dim(x = vector2, KernelsOrNbins = xA)
  y2 = data.frame(dens2)
  dfPDF <- merge(y1, y2, by = "Kernels")
  EucDist <- dist(rbind(dfPDF$SDH.x, dfPDF$SDH.y))
  return(EucDist)
}
