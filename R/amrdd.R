#Calculates the mean relative difference for two probability density distributions
#expressed as Smooth Density Histogram
#' @importFrom SmoothDensHist1dim SmoothDensHist1dim
amrdd <- function(vector1, vector2) {
  densA <- SmoothDensHist1dim(c(vector1,vector2))
  xA = densA$Kernels
  dens1 <- SmoothDensHist1dim(x = vector1, KernelsOrNbins = xA)
  y1 = data.frame(dens1)
  dens2 <- SmoothDensHist1dim(x = vector2, KernelsOrNbins = xA)
  y2 = data.frame(dens2)
  dfPDF <- merge(y1,y2, by = "Kernels")
  dfPDF$Diff <- dfPDF$SDH.x - dfPDF$SDH.y
    dfPDF$Mean <- (dfPDF$SDH.x + dfPDF$SDH.y) / 2
  dfPDF$Mean[dfPDF$Mean == 0] <- 1
  dfPDF$RelDiff <- dfPDF$Diff / dfPDF$Mean
  RelDiff <- abs(mean(dfPDF$RelDiff, na.rm = TRUE))
  return(RelDiff)
}
