#' Calculate a Smooth Density Histogram for a Single Variable
#'
#' This function calculates a smooth density histogram for a single variable.
#'
#' @param x A numeric vector representing the data.
#' @param KernelsOrNbins Either the number of bins or a vector of kernel positions.
#' @param SDHinPercent A logical value indicating whether the Smooth Density Histogram
#'   should be returned as a percentage (TRUE) or not (FALSE).
#' @param lambda The smoothing parameter for the Smooth Density Histogram.
#'
#' @return A list with two elements: `Kernels` (the kernel positions) and `SDH` (the
#'   Smooth Density Histogram values).
#'
#' @importFrom pracma histc
#' @importFrom caTools trapz
#'
SmoothDensHist1dim <- function(x, KernelsOrNbins = NULL, SDHinPercent = FALSE, lambda = 20) {
  # Check if the input vector is empty
  if (length(x) == 0) {
    warning("SmoothDensHist1dim: Size of x is zero.", call. = FALSE)
    if (is.null(KernelsOrNbins)) {
      Kernels <- 1
    } else {
      Kernels <- KernelsOrNbins
    }
    return(list(Kernels = Kernels, SDH = Kernels * 0))
  }
  
  # Load required packages
  requireNamespace("pracma")
  requireNamespace("caTools")
  
  # Define a helper function for smoothing the histogram
  smooth1D <- function(Y, lambda) {
    if (is.vector(Y)) {
      Y <- as.matrix(Y)
    }
    Y[is.na(Y)] <- 0
    dd <- dim(Y)
    m <- dd[1]
    E <- diag(m)
    D1 <- (diff(E, 1))
    D2 <- (diff(D1, 1))
    P <- lambda^2 * (t(D2) %*% D2) + 2 * lambda * (t(D1) %*% D1)
    Z <- solve((E + P), Y)
    return(Z)
  }
  
  # Set default values for missing parameters
  if (missing(lambda)) {
    lambda <- 20
  }
  if (missing(SDHinPercent)) {
    SDHinPercent <- FALSE
  }
  if (is.null(KernelsOrNbins)) {
    KernelsOrNbins <- 200
  }
  if (length(KernelsOrNbins) < 1) {
    KernelsOrNbins <- 200
  }
  
  # Remove non-finite values from the input vector
  x <- x[is.finite(x)]
  n <- length(x)
  minx <- min(x, na.rm = TRUE)
  maxx <- max(x, na.rm = TRUE)
  
  # Determine the kernel positions and bin edges
  if (length(KernelsOrNbins) == 1) {
    nbins <- KernelsOrNbins
    edges1 <- seq(from = minx, to = maxx, length.out = (nbins + 1))
    end <- length(edges1)
    Kernels <- c(edges1[1:(end - 1)] + 0.5 * diff(edges1))
    InInd <- c()
  } else {
    Kernels <- c(KernelsOrNbins)
    InInd <- which((Kernels >= minx) & (Kernels <= maxx))
  }
  
  # Calculate the Smooth Density Histogram
  if (length(InInd) == 0) {
    SDH <- Kernels * 0
    nbins <- 1
  } else {
    DataInd <- which((x >= min(Kernels[InInd])) & (x <= max(Kernels[InInd])))
    if (length(DataInd) < 2) {
      if (length(DataInd) == 0) {
        SDH <- Kernels * 0
        nbins <- 1
      } else {
        SDH <- Kernels * 0
        SDH[InInd] <- 1
        nbins <- 1
      }
    } else {
      x <- x[DataInd]
      edges1 <- Kernels[InInd]
      nbins <- length(edges1)
      edges1 <- c(-Inf, edges1[2:(end - 1)], Inf)
      V <- pracma::histc(x, edges1)
      dummy <- V$cnt
      H <- dummy / n
      SDH <- smooth1D(H, nbins / lambda)
      SDH <- as.vector(SDH)
    }
    
    if (length(DataInd) > 1) {
      sdh <- SDH
      SDH <- Kernels * 0
      SDH[InInd] <- sdh
    }
  }
  
  # Normalize the Smooth Density Histogram
  if (SDHinPercent) {
    SDH <- SDH / max(SDH)
  } else {
    if (sum(SDH) == 0) {
      Area <- 0
    } else {
      Area <- caTools::trapz(Kernels, SDH)
    }
    if (Area < 1e-10) {
      SDH <- rep(0, length(Kernels))
    } else {
      SDH <- SDH / Area
    }
  }
  
  return(list(Kernels = Kernels, SDH = SDH))
}
