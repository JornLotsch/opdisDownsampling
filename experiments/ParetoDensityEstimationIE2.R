#' Calculate Optimal Number of Bins Using Izenman-Einstein Rule
#'
#' This function calculates the optimal number of bins for histogram construction
#' using a variant of the Izenman-Einstein rule, which balances bias and variance
#' in density estimation.
#'
#' @param Data A numeric vector, matrix, or data frame containing the data.
#'   If matrix, each column is treated separately. Non-finite values are handled
#'   automatically.
#'
#' @return An integer representing the optimal number of bins. Returns 0 if no
#'   valid data is available, and at least 10 bins for valid data.
#'
#' @details The function implements a robust bin selection method that:
#' \itemize{
#'   \item Uses Silverman's rule of thumb with interquartile range adjustment
#'   \item Handles both small (< 5000) and large datasets with appropriate quantile methods
#'   \item Provides a minimum of 10 bins for meaningful histograms
#'   \item Uses different quantile estimation methods based on data size for efficiency
#' }
#'
#' @references
#' Izenman, A.J. and Sommer, C.J. (1988). Philatelic mixtures and multimodal densities.
#' Journal of the American Statistical Association, 83, 941-953.
#'
#' @importFrom stats sd quantile
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' data <- rnorm(1000)
#' nbins <- optimal_no_bins_ie(data)
#'
#' # With matrix input
#' mat_data <- matrix(rnorm(2000), ncol = 2)
#' nbins <- optimal_no_bins_ie(mat_data)
#' }
#'
#' @export
optimal_no_bins_ie <- function(Data) {
  # Determine data size and handle different input types
  if (is.matrix(Data)) {
    nData <- colSums(!is.nan(Data))
  } else if (is.vector(Data)) {
    nData <- sum(!is.nan(Data))
  } else {
    nData <- 0
  }

  # Handle insufficient data
  if (nData < 1) {
    optNrOfBins <- 0
  } else {
    # Calculate standard deviation
    sigma <- stats::sd(Data, na.rm = TRUE)

    # Use appropriate quantile method based on data size
    if (nData < 5000) {
      p <- c_quantile(na.omit(Data), c(0.25, 0.75))
    } else {
      p <- stats::quantile(Data, c(0.25, 0.75), type = 8, na.rm = TRUE)
    }

    # Calculate robust scale estimate using interquartile range
    interquartilRange <- p[2] - p[1]
    sigmaSir <- min(sigma, interquartilRange / 1.349)

    # Apply Izenman-Einstein rule
    optBinWidth <- 3.49 * sigmaSir / (nData)^(1 / 3)

    # Convert bin width to number of bins
    if (optBinWidth > 0) {
      dataRange <- max(Data, na.rm = TRUE) - min(Data, na.rm = TRUE)
      optNrOfBins <- max(ceiling(dataRange / optBinWidth), 10)
    } else {
      optNrOfBins <- 10
    }
  }

  return(optNrOfBins)
}

#' Calculate Pareto Radius for Density Estimation
#'
#' This function computes the Pareto radius used in Pareto density estimation,
#' which represents the optimal bandwidth for kernel-based density estimation
#' in high-dimensional spaces.
#'
#' @param Data A numeric vector containing the data points for radius calculation.
#' @param maximumNrSamples An integer specifying the maximum number of samples
#'   to use in distance calculations (default: 10000). For large datasets,
#'   sampling reduces computational cost.
#'
#' @return A numeric value representing the Pareto radius. The radius is
#'   automatically adjusted for large datasets to maintain estimation quality.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item Validates input data for NaN, NA, and infinite values
#'   \item Samples data if the dataset exceeds maximumNrSamples
#'   \item Calculates pairwise Euclidean distances
#'   \item Selects the 18th percentile as the initial radius
#'   \item Adjusts radius for large datasets (> 1024 observations)
#' }
#'
#' @section Warning:
#' The function will issue warnings for problematic data and will stop execution
#' if the radius cannot be calculated due to NaN, NA, or infinite results.
#'
#' @importFrom stats runif quantile
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' data <- rnorm(500)
#' radius <- pareto_radius_ie(data)
#'
#' # Large dataset with sampling
#' large_data <- rnorm(50000)
#' radius <- pareto_radius_ie(large_data, maximumNrSamples = 5000)
#' }
#'
#' @export
pareto_radius_ie <- function(Data, maximumNrSamples = 10000) {
  requireNamespace("parallelDist")

  # Validate input data quality
  ntemp <- sum(is.nan(Data))
  if (ntemp > 0) {
    warning("Data contains NaN values, Pareto Radius may not be calculated correctly.", call. = FALSE)
  }

  ntemp2 <- sum(is.na(Data))
  if (ntemp2 > ntemp) {
    warning("Data contains NA values, Pareto Radius may not be calculated correctly.", call. = FALSE)
  }

  ntemp3 <- sum(is.infinite(Data))
  if (ntemp3 > 0) {
    warning("Data contains infinite values, Pareto Radius may not be calculated correctly.", call. = FALSE)
  }

  nData <- length(Data)

  # Sample data if necessary to reduce computational cost
  if (maximumNrSamples >= nData) {
    sampleData <- Data
  } else {
    sampleInd <- ceiling(stats::runif(maximumNrSamples, min = 0, max = nData))
    sampleData <- Data[sampleInd]
  }

  # Calculate pairwise distances efficiently
  distvec <- as.vector(parallelDist::parallelDist(
    as.matrix(sampleData),
    method = "euclidean",
    upper = FALSE,
    diag = FALSE
  ))

  # Select Pareto radius using 18th percentile
  if (nData < 5000) {
    paretoRadius <- c_quantile(na.omit(distvec), 18 / 100)
  } else {
    paretoRadius <- stats::quantile(distvec, 18 / 100, type = 8, na.rm = TRUE)
  }

  # Handle zero radius case
  if (paretoRadius == 0) {
    if (nData < 5000) {
      pzt <- c_quantile(na.omit(distvec), probs = c(1:100) / 100)
    } else {
      pzt <- stats::quantile(distvec, probs = c(1:100) / 100, type = 8, na.rm = TRUE)
    }
    paretoRadius <- min(pzt[pzt > 0], na.rm = TRUE)
  }

  # Validate radius calculation results
  if (is.nan(paretoRadius)) {
    stop("Pareto Radius could not be calculated (NaN result). Check input data quality.", call. = FALSE)
  }
  if (is.na(paretoRadius)) {
    stop("Pareto Radius could not be calculated (NA result). Check input data quality.", call. = FALSE)
  }
  if (!is.finite(paretoRadius)) {
    stop("Pareto Radius could not be calculated (infinite result). Check input data quality.", call. = FALSE)
  }

  # Adjust radius for large datasets to maintain estimation quality
  if (nData > 1024) {
    paretoRadius <- paretoRadius * 4 / (nData^0.2)
  }

  return(paretoRadius)
}

#' Pareto Density Estimation for One-Dimensional Data
#'
#' This function performs Pareto density estimation on univariate data using
#' kernel-based methods with automatic boundary correction and adaptive
#' kernel placement.
#'
#' @param Data A numeric vector containing the data for density estimation.
#'   Non-finite values are automatically removed.
#' @param paretoRadius A numeric value specifying the Pareto radius (bandwidth).
#'   If missing or invalid, it will be calculated automatically using
#'   \code{pareto_radius_ie()}.
#' @param kernels A numeric vector specifying kernel positions, or NULL for
#'   automatic kernel placement based on data range and optimal binning.
#' @param MinAnzKernels An integer specifying the minimum number of kernels
#'   to use (default: 100). Used when kernels are generated automatically.
#'
#' @return A list with three elements:
#'   \item{kernels}{Numeric vector of kernel positions}
#'   \item{paretoDensity}{Numeric vector of density estimates at kernel positions}
#'   \item{paretoRadius}{The Pareto radius used in the estimation}
#'
#' @details The function implements several key features:
#' \itemize{
#'   \item Automatic data validation and cleaning
#'   \item Special handling for degenerate cases (< 3 unique values)
#'   \item Boundary correction using data reflection
#'   \item Automatic kernel placement using optimal binning
#'   \item Density normalization to ensure proper probability distribution
#' }
#'
#' @section Special Cases:
#' For datasets with fewer than 3 unique values, the function assumes
#' Dirac delta distributions and creates appropriate spike densities.
#'
#' @importFrom caTools trapz
#'
#' @examples
#' \dontrun{
#' # Basic density estimation
#' data <- rnorm(1000)
#' result <- pareto_density_estimation_ie(data)
#' plot(result$kernels, result$paretoDensity, type = "l")
#'
#' # Custom radius and kernels
#' custom_kernels <- seq(-3, 3, length.out = 50)
#' result <- pareto_density_estimation_ie(data, paretoRadius = 0.5,
#'                                        kernels = custom_kernels)
#' }
#'
#' @export
pareto_density_estimation_ie <- function(Data, paretoRadius, kernels = NULL, MinAnzKernels = 100) {
  requireNamespace("caTools")

  # Validate and convert data to numeric vector
  if (!is.vector(Data)) {
    Data <- as.vector(Data)
    warning("Data converted to vector. Ensure data is univariate for proper density estimation.", call. = FALSE)
  }

  if (!is.numeric(Data)) {
    Data <- as.numeric(Data)
    warning("Data converted to numeric. Check for potential data loss during conversion.", call. = FALSE)
  }

  # Handle non-finite values
  if (length(Data) != sum(is.finite(Data))) {
    message("Non-finite values detected and will be removed from density estimation.")
  }
  Data <- Data[is.finite(Data)]
  values <- unique(Data)

  # Handle degenerate cases with few unique values
  if (length(values) > 2 && length(values) < 5) {
    warning("Fewer than 5 unique values detected. Density estimation may be unreliable.", call. = FALSE)
  }

  if (length(values) < 3) {
    warning("Fewer than 3 unique values detected. Assuming Dirac delta distribution(s).",
            call. = FALSE)

    # Create delta function approximation for first value
    if (values[1] != 0) {
      kernels <- seq(from = values[1] * 0.9, to = values[1] * 1.1, by = values[1] * 0.0001)
    } else {
      kernels <- seq(from = values[1] - 0.1, to = values[1] + 0.1, by = 0.0001)
    }

    paretoDensity <- rep(0, length(kernels))
    paretoDensity[kernels == values[1]] <- 1

    # Handle second value if present
    if (length(values) == 2) {
      if (values[2] != 0) {
        kernels2 <- seq(from = values[2] * 0.9, to = values[2] * 1.1, by = values[2] * 0.0001)
      } else {
        kernels2 <- seq(from = values[2] - 0.1, to = values[2] + 0.1, by = 0.0001)
      }

      paretoDensity2 <- rep(0, length(kernels2))
      paretoDensity2[kernels2 == values[2]] <- 1
      paretoDensity <- c(paretoDensity, paretoDensity2)
      kernels <- c(kernels, kernels2)
    }

    return(list(kernels = kernels, paretoDensity = paretoDensity, paretoRadius = 0))
  }

  # Warn about potential issues with small datasets
  if (length(Data) < 10) {
    warning("Fewer than 10 data points provided. Radius calculation may be unreliable.", call. = FALSE)
  }

  # Calculate Pareto radius if not provided or invalid
  if (missing(paretoRadius) || is.null(paretoRadius) || is.na(paretoRadius) ||
    paretoRadius == 0 || length(paretoRadius) == 0) {
    paretoRadius <- pareto_radius_ie(Data)
  }

  # Determine data range
  minData <- min(Data, na.rm = TRUE)
  maxData <- max(Data, na.rm = TRUE)

  # Generate kernels automatically if not provided
  if (length(kernels) <= 1) {
    if (length(kernels) == 0 || (length(kernels) == 1 && kernels == 0)) {
      # Use optimal binning to determine kernel positions
      nBins <- optimal_no_bins_ie(Data)
      nBins <- max(MinAnzKernels, nBins)

      # Handle excessive bin numbers
      if (nBins > 100) {
        if (nBins > 1E4) {
          nBins <- 1E4
          warning("Excessive number of bins estimated. Consider data transformation or sampling.", call. = FALSE)
        } else {
          nBins <- nBins * 3 + 1
        }
      }

      # Generate kernel positions using pretty breaks
      breaks <- pretty(c(minData, maxData), n = nBins, min.n = 1)
      nB <- length(breaks)
      mids <- 0.5 * (breaks[-1L] + breaks[-nB])
      kernels <- mids
    }
  }

  nKernels <- length(kernels)

  # Apply boundary correction using data reflection
  lowBInd <- (Data < (minData + paretoRadius))
  lowR <- as.matrix(2 * minData - Data[lowBInd], ncol = 1)
  upBInd <- (Data > (maxData - paretoRadius))
  upR <- as.matrix(2 * maxData - Data[upBInd], ncol = 1)
  DataPlus <- as.matrix(c(Data, lowR, upR), 1)

  # Calculate density estimates at each kernel position
  paretoDensity <- rep(0, nKernels)
  for (i in 1:nKernels) {
    lb <- kernels[i] - paretoRadius
    ub <- kernels[i] + paretoRadius
    isInParetoSphere <- (DataPlus >= lb) & (DataPlus <= ub)
    paretoDensity[i] <- sum(isInParetoSphere)
  }

  # Normalize density to ensure proper probability distribution
  area <- caTools::trapz(kernels, paretoDensity)
  if (area < 1e-10 || is.na(area)) {
    paretoDensity <- rep(0, nKernels)
  } else {
    paretoDensity <- paretoDensity / area
  }

  return(list(kernels = kernels, paretoDensity = paretoDensity, paretoRadius = paretoRadius))
}