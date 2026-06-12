#' Perform Class-Proportional Downsampling
#'
#' This function samples a subset of data based on the similarity of its probability
#' distribution to that of the original data.
#'
#' @param Data A data frame containing the data to be downsampled.
#' @param Cls A vector of class labels for the data. If not provided, all instances
#'   will be assigned to a single class.
#' @param Size A numeric value specifying the desired size of the downsampled dataset.
#'   If 0 < Size < 1, it is treated as a proportion of the original dataset.
#'   If Size >= 1, it is treated as the absolute number of instances to retain.
#' @param Seed Seed value. Can be an integer value, or one of:
#'   \itemize{
#'     \item \code{"auto"}: Uses seed recovery to match the current RNG state.
#'     \item \code{"simple"} (default): Generates and reports a seed using the
#'       current RNG state.
#'     \item Integer: Uses the supplied seed value for exact reproducibility.
#'   }
#'   For systematic testing or fully reproducible analyses, use integer values.
#' @param nTrials The number of trials to perform when sampling the data.
#' @param TestStat A character string specifying the statistical test to be used for
#'   comparing the distributions.
#' @param MaxCores The maximum number of cores to use for parallel processing.
#' @param PCAimportance A logical value indicating whether to use PCA to identify
#'   relevant variables.
#' @param JobSize Number of seeds to process in each chunk for memory optimization.
#'   If \code{0}, no chunking is applied. If \code{NULL}, an automatic chunk size
#'   is determined based on data size, nTrials, and available memory.
#' @param verbose Logical, whether to print chunk size diagnostics.
#'
#' @return A list with the following elements:
#'   - `ReducedData`: The downsampled dataset.
#'   - `RemovedData`: The removed data from the original dataset.
#'   - `ReducedInstances`: The row names of the downsampled dataset.
#'   - `RemovedInstances`: The row names of the removed dataset.
#'
#' @useDynLib opdisDownsampling, .registration = TRUE
#' @importFrom methods hasArg
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom pbmcapply pbmclapply
#' @import foreach
#' @export
opdisDownsampling <- function(Data, Cls, Size, Seed = "simple", nTrials = 1000, TestStat = "ad",
                              MaxCores = getOption("mc.cores", 2L), PCAimportance = FALSE,
                              JobSize = 0, verbose = FALSE) {
  # Create empty data frame
  dfx <- data.frame(Data)
  dfxempty <- dfx[0, ]

  # Check if correct input is provided and library can be run
  if (!all(sapply(dfx, is.numeric))) {
    stop("opdisDownsampling: Only numeric data allowed. Nothing to downsample.")
  }

  if (!is.numeric(Size) || length(Size) != 1 || is.na(Size) || !is.finite(Size)) {
    stop("opdisDownsampling: Size must be a single finite numeric value.")
  }

  if (!is.numeric(nTrials) || length(nTrials) != 1 || is.na(nTrials) || !is.finite(nTrials) || nTrials < 1) {
    stop("opdisDownsampling: nTrials must be a positive integer.")
  }
  nTrials <- as.integer(nTrials)

  # Check for variables with excessive NA values that might cause issues with small subsamples
  data_cols <- names(dfx)
  for (col in data_cols) {
    na_proportion <- sum(is.na(dfx[[col]])) / nrow(dfx)

    # Warn if a variable has all NAs
    if (na_proportion == 1) {
      warning(sprintf(
        "opdisDownsampling: Variable '%s' contains only NA values and will be excluded from distribution comparisons.",
        col
      ), call. = FALSE)
    }
    # Warn if a variable has many NAs and subsample is small (before Size conversion)
    else if (na_proportion > 0.5) {
      # Check if Size is a proportion or absolute
      target_size <- if (Size > 0 && Size < 1) Size * nrow(dfx) else Size
      if (target_size < nrow(dfx) * 0.1) {
        warning(sprintf(
          "opdisDownsampling: Variable '%s' has %.1f%% NA values. With small subsample size (~%d), there is risk of drawing only NAs.",
          col, na_proportion * 100, round(target_size)
        ), call. = FALSE)
      }
    }
  }

  # Handle class labels
  Clsarg_missing <- missing(Cls)
  if (!Clsarg_missing) {
    if (length(Cls) != nrow(dfx)) {
      stop("opdisDownsampling: Unequal number of cases and class memberships.")
    }
  } else {
    Cls <- rep(1, nrow(dfx))
  }
  dfx$Cls <- Cls

  if (!is.numeric(Size) || length(Size) != 1 || is.na(Size) || !is.finite(Size)) {
    stop("opdisDownsampling: Size must be a single finite numeric value.")
  }

  if (!is.numeric(nTrials) || length(nTrials) != 1 || is.na(nTrials) || !is.finite(nTrials) || nTrials < 1) {
    stop("opdisDownsampling: nTrials must be a positive integer.")
  }
  nTrials <- as.integer(nTrials)

  # Handle size parameter: convert proportion to absolute count if needed
  original_size <- Size # Keep original for potential error messages
  if (Size > 0 && Size < 1) {
    # Convert proportion to absolute count
    Size <- round(Size * nrow(dfx))
    if (verbose) {
      message(sprintf("Size converted from proportion %.3f to absolute count: %d", original_size, Size))
    }
  } else if (Size >= 1) {
    # Use as absolute count (ensure it's integer)
    Size <- as.integer(Size)
  }

  # Validate the final size
  if (Size >= nrow(dfx)) {
    warning(sprintf(
      "opdisDownsampling: Size (%d) >= number of rows (%d). Nothing to downsample.",
      Size, nrow(dfx)
    ), call. = FALSE)

    ReducedData <- dfx
    RemovedData <- dfxempty

    if (Clsarg_missing) {
      ReducedData <- ReducedData[1:(ncol(ReducedData) - 1)]
      RemovedData <- RemovedData[1:(ncol(RemovedData) - 1)]
    }

    return(list(
      ReducedData = ReducedData,
      RemovedData = RemovedData,
      ReducedInstances = rownames(ReducedData),
      RemovedInstances = rownames(RemovedData)
    ))
  }

  if (Size <= 0) {
    warning(sprintf("opdisDownsampling: Size (%d) <= 0. All data will be removed.", Size), call. = FALSE)

    ReducedData <- dfxempty
    RemovedData <- dfx

    if (Clsarg_missing) {
      ReducedData <- ReducedData[1:(ncol(ReducedData) - 1)]
      RemovedData <- RemovedData[1:(ncol(RemovedData) - 1)]
    }

    return(list(
      ReducedData = ReducedData,
      RemovedData = RemovedData,
      ReducedInstances = rownames(ReducedData),
      RemovedInstances = rownames(RemovedData)
    ))
  }

  # Handle test statistic
  TestStats <- c("ad", "kuiper", "cvm", "wass", "dts", "ks", "kld", "amrdd", "euc", "nent")
  if (!(TestStat %in% TestStats)) {
    warning(paste0(
      "opdisDownsampling: Possible TestStat = ", paste(TestStats, collapse = ", "),
      ". TestStat set to default = 'ad'."
    ), call. = FALSE)
    TestStat <- "ad"
  }

  # Initialize environment
  # Seed handling with three clear options
  if (is.numeric(Seed)) {
    # Option 1: Use provided integer seed directly
    Seed <- as.integer(Seed)
  } else if (is.character(Seed)) {
    Seed <- switch(Seed,
      "auto" = as.integer(get_seed()), # Complex seed recovery
      "simple" = { # Generate and report a seed from the current RNG state
        temp_seed <- sample(1:100000, 1)
        warning(paste0("opdisDownsampling: Seed set at ", temp_seed, "."), call. = FALSE)
        temp_seed
      },
      stop("Invalid Seed input. Use 'auto', 'simple', or an integer.")
    )
  } else {
    # Fallback for backward compatibility
    Seed <- as.integer(get_seed())
  }

  list.of.seeds <- 1:nTrials + Seed - 1

  # Number of cores used
  nProc <- determine_n_cores(MaxCores)

  # Determine chunk size
  if (is.null(JobSize)) {
    # NULL enables automatic memory-aware chunk-size calculation
    JobSize <- calculate_optimal_chunk_size(
      n_rows = nrow(dfx),
      n_cols = ncol(dfx) - 1,
      nTrials = nTrials,
      nProc = nProc
    )

    print_chunk_diagnostics(nrow(dfx), ncol(dfx) - 1, nTrials, JobSize, verbose)
  } else if (JobSize <= 0) {
    # 0 or negative means no chunking: process all trials in one batch
    JobSize <- nTrials

    if (verbose) {
      data_size_mb <- (nrow(dfx) * (ncol(dfx) - 1) * 8) / (1024^2)
      message(sprintf("Chunk size diagnostics:"))
      message(sprintf("  Data: %d rows x %d cols (%.1f MB)", nrow(dfx), ncol(dfx) - 1, data_size_mb))
      message(sprintf(
        "  Trials: %d, Chunk size: %d (no chunking - single batch)",
        nTrials, JobSize
      ))
    }
  } else {
    # Positive JobSize means user-defined chunk size
    JobSize <- as.integer(JobSize)
    print_chunk_diagnostics(nrow(dfx), ncol(dfx) - 1, nTrials, JobSize, verbose)
  }

  # Perform sampling and analyze picked data subsets
  ReducedDiag <- sample_and_analyze(
    DataAndClasses = dfx,
    TestStat = TestStat,
    Size = Size,
    list.of.seeds = list.of.seeds,
    PCAimportance = PCAimportance,
    nProc = nProc,
    JobSize = JobSize
  )

  # Validate and process trial results
  validate_reduced_diag(ReducedDiag)

  # Memory-efficient matrix construction with enhanced error handling
  n_trials <- length(ReducedDiag)

  # Enhanced validation of data structure
  if (!is.list(ReducedDiag[[1]]) ||
    is.null(ReducedDiag[[1]][[1]]) ||
    length(ReducedDiag[[1]][[1]]) == 0) {
    stop("opdisDownsampling: Invalid structure in trial results.")
  }

  n_vars <- length(ReducedDiag[[1]][[1]])
  var_names <- names(ReducedDiag[[1]][[1]])

  # Pre-allocate matrices with error handling
  tryCatch(
    {
      AD_reduced_statMat <- matrix(NA_real_,
        nrow = n_trials, ncol = n_vars,
        dimnames = list(NULL, var_names)
      )
    },
    error = function(e) {
      stop(sprintf("opdisDownsampling: Failed to allocate matrices: %s", e$message))
    }
  )

  # Fill matrices and validate
  for (i in seq_len(n_trials)) {
    if (!is.null(ReducedDiag[[i]])) {
      tryCatch(
        {
          # Validate data structure for each trial
          if (is.list(ReducedDiag[[i]]) &&
            all(c("ADv_reduced", "ADv_removed") %in% names(ReducedDiag[[i]]))) {
            AD_reduced_statMat[i, ] <- ReducedDiag[[i]][["ADv_reduced"]]
          } else {
            warning(sprintf("opdisDownsampling: Invalid data structure in trial %d, skipping.", i),
              call. = FALSE
            )
          }
        },
        error = function(e) {
          warning(sprintf("opdisDownsampling: Error processing trial %d: %s", i, e$message),
            call. = FALSE
          )
        }
      )
    }
  }

  # Clean up the large ReducedDiag object early
  rm(ReducedDiag)
  gc()

  # Validate matrices before proceeding
  validate_matrices(
    list(AD_reduced_statMat),
    c("AD_reduced_statMat")
  )

  # Find best subsample using refactored selection logic with error handling
  BestTrial <- tryCatch(
    {
      select_best_trial_one_matrix(AD_reduced_statMat)
    },
    error = function(e) {
      warning(sprintf("opdisDownsampling: Error in trial selection: %s. Using first trial.", e$message),
        call. = FALSE
      )
      1
    }
  )

  # Validate BestTrial result
  if (is.na(BestTrial) || BestTrial < 1 || BestTrial > n_trials) {
    warning(sprintf(
      "opdisDownsampling: Invalid best trial index (%s). Using first trial.",
      as.character(BestTrial)
    ), call. = FALSE)
    BestTrial <- 1
  }

  # Clean up stat matrices - moved after BestTrial validation
  rm(AD_reduced_statMat)
  gc()

  # Get the best subsample
  df_reduced_final <- MakeReducedDataMat(DataAndClasses = dfx, Size = Size, Seed = list.of.seeds[BestTrial])
  ReducedData <- df_reduced_final$ReducedDataList
  RemovedData <- df_reduced_final$RemovedDataList

  if (Clsarg_missing) {
    ReducedData <- ReducedData[1:(ncol(ReducedData) - 1)]
    RemovedData <- RemovedData[1:(ncol(RemovedData) - 1)]
  }

  return(list(
    ReducedData = ReducedData,
    RemovedData = RemovedData,
    ReducedInstances = rownames(ReducedData),
    RemovedInstances = rownames(RemovedData)
  ))
}
