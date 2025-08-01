#' Calculate a Smooth Density Histogram for a Single Variable
#'
#' This function computes a smooth density histogram using regularized smoothing
#' for a given dataset. It handles both user-specified kernel positions and
#' automatic bin generation, with optional percentage-based output.
#'
#' @param x A numeric vector representing the input data.
#' @param KernelsOrNbins Either the number of bins (integer) or a vector of
#'   kernel positions (numeric vector). If NULL, defaults to 200 bins.
#' @param SDHinPercent A logical value indicating whether the output should be
#'   normalized as percentages (TRUE) or as a probability density (FALSE).
#' @param lambda The regularization parameter controlling smoothness. Higher
#'   values produce smoother results. Default is 20.
#'
#' @return A list with two elements:
#'   \item{Kernels}{Numeric vector of kernel positions}
#'   \item{SDH}{Numeric vector of smooth density histogram values}
#'
#' @importFrom pracma histc linspace
#' @importFrom caTools trapz
#'
SmoothDensHist1dim <- function(x, KernelsOrNbins = NULL, SDHinPercent = FALSE, lambda = 20) {
  # Load required packages
  requireNamespace("pracma")
  requireNamespace("caTools")

  # Helper function for 1D regularized smoothing using finite differences
  smooth1D <- function(Y, lambda) {
    if (is.vector(Y)) {
      Y <- as.matrix(Y)
    }
    Y[is.na(Y)] <- 0
    dd <- dim(Y)
    m <- dd[1]

    # Create identity matrix and difference operators
    E <- diag(m)
    D1 <- diff(E, 1) # First difference operator
    D2 <- diff(D1, 1) # Second difference operator

    # Regularization matrix combining first and second order penalties
    P <- lambda ^ 2 * (t(D2) %*% D2) + 2 * lambda * (t(D1) %*% D1)

    # Solve regularized system
    Z <- solve((E + P), Y)
    return(Z)
  }

  # Handle empty input data
  if (length(x) == 0) {
    warning("SmoothDensHist1D: Size of x is zero.", call. = FALSE)
    if (is.null(KernelsOrNbins)) {
      Kernels <- 1
    } else {
      Kernels <- KernelsOrNbins
    }
    SDH <- Kernels * 0
    return(list(Kernels = Kernels, SDH = SDH))
  }

  # Set default parameters if missing
  if (missing(lambda)) {
    lambda <- 20
  }
  if (missing(SDHinPercent)) {
    SDHinPercent <- FALSE
  }

  # Set default number of bins
  if (is.null(KernelsOrNbins)) {
    KernelsOrNbins <- 200
  }
  if (length(KernelsOrNbins) < 1) {
    KernelsOrNbins <- 200
  }

  # Remove non-finite values from input data
  x <- x[is.finite(x)]
  n <- length(x)
  minx <- min(x, na.rm = TRUE)
  maxx <- max(x, na.rm = TRUE)

  # Generate kernels based on input type
  if (length(KernelsOrNbins) == 1) {
    # Number of bins specified - generate equally spaced kernels
    nbins <- KernelsOrNbins
    edges1 <- pracma::linspace(minx, maxx, nbins + 1)
    end <- length(edges1)
    # Use bin centers as kernel positions
    Kernels <- edges1[1:(end - 1)] + 0.5 * diff(edges1)
    InInd <- c()
  } else {
    # Kernel positions specified directly
    Kernels <- c(KernelsOrNbins)
  }

  # Find kernels within data range
  InInd <- which((Kernels >= minx) & (Kernels <= maxx))

  if (length(InInd) == 0) {
    # No kernels in data range
    SDH <- Kernels * 0
    nbins <- 1
  } else {
    # Find data points within kernel range
    DataInd <- which((x >= min(Kernels[InInd])) & (x <= max(Kernels[InInd])))

    if (length(DataInd) < 2) {
      # Insufficient data points
      if (length(DataInd) == 0) {
        SDH <- Kernels * 0
        nbins <- 1
      } else {
        SDH <- Kernels * 0
        SDH[InInd] <- 1
        nbins <- 1
      }
    } else {
      # Process data with sufficient points
      x <- x[DataInd]
      edges1 <- sort(Kernels[InInd], decreasing = FALSE)
      nbins <- length(edges1)

      # Prepare histogram edges for pracma::histc
      ende <- length(edges1)
      if (length(edges1) > 3) {
        edges1 <- c(-Inf, edges1[2:(ende - 1)], Inf)
      } else if (length(edges1) == 3) {
        edges1 <- c(-Inf, edges1[2], Inf)
      } else {
        warning("SmoothDensHist1D: Two or fewer kernels in data range. Density estimation may be inaccurate.")
        edges1 <- c(-Inf, edges1, Inf)
      }

      # Calculate histogram
      V <- pracma::histc(x, edges = edges1)
      dummy <- V$cnt
      bin <- V$bin

      # Normalize to get density estimate
      H <- dummy / n

      # Apply smoothing
      SDH <- smooth1D(H, nbins / lambda)
      SDH <- as.vector(SDH)
    }

    # Map smoothed values back to full kernel range if needed
    if (length(InInd) > 1) {
      sdh <- SDH
      SDH <- Kernels * 0
      SDH[InInd] <- sdh
    }
  }

  # Normalize output based on requested format
  if (!SDHinPercent) {
    # Normalize as probability density (area under curve = 1)
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
  } else {
    # Convert to percentage scale
    SDH <- SDH / nanmax(SDH)
  }

  return(list(Kernels = Kernels, SDH = SDH))
}