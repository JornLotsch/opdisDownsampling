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
#' @param Seed An integer value to be used as the random seed for reproducibility.
#' @param nTrials The number of trials to perform when sampling the data.
#' @param TestStat A character string specifying the statistical test to be used for
#'   comparing the distributions.
#' @param MaxCores The maximum number of cores to use for parallel processing.
#' @param PCAimportance A logical value indicating whether to use PCA to identify
#'   relevant variables.
#' @param CheckRemoved A logical value indicating whether to also optimize the removed part
#'   of the data for distribution equality with the original.
#' @param CheckThreefold A logical value indicating whether to also optimize the reduced part
#'   of the data for distribution equality with the removed part. Ignored when CheckRemoved is FALSE.
#' @param OptimizeBetween A logical value indicating whether to optimize the reduced part
#'   of the data for distribution equality with the removed part. If set, all other comparisons are not performed.
#' @param JobSize Number of seeds to process in each chunk for memory optimization.
#'   If NULL, automatically determined based on data size, nTrials, and available memory.
#' @param verbose Logical, whether to print chunk size diagnostics.
#' @param NonNoiseSelection A logical value indicating whether to use non-uniform
#'   distribution tests to identify relevant variables.
#' @param UniformTestStat A character string specifying the statistical test to be used for
#'   non-uniform variable selection. Available options are: "ks" (Kolmogorov-Smirnov),
#'   "ad" (Anderson-Darling), "kuiper", "cvm" (Cramér-von Mises), "wass" (Wasserstein),
#'   "dts" (Distributional Transform Statistic), "kld" (Kullback-Leibler divergence),
#'   "amrdd" (Average Mean Root of Distributional Differences), and "euc" (Euclidean distance).
#'   Only used when NonNoiseSelection = TRUE (default: "ks").
#' @param UniformThreshold Threshold value for non-uniform variable selection (default: 0.1).
#' @param WorstSample A logical value for testing purpose reversing the split ranking to obtain
#'   the least similar subsample (default: FALSE).
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
opdisDownsampling <- function(Data, Cls, Size, Seed, nTrials = 1000, TestStat = "ad",
                              MaxCores = getOption("mc.cores", 2L), PCAimportance = FALSE,
                              CheckRemoved = FALSE, CheckThreefold = FALSE, OptimizeBetween = FALSE,
                              JobSize = 0, verbose = FALSE, NonNoiseSelection = FALSE,
                              UniformTestStat = "ks", UniformThreshold = 0.05, WorstSample = FALSE) {

  # Set CheckThreefold to FALSE when CheckRemoved is FALSE
  if (!CheckRemoved) CheckThreefold <- FALSE

  # Set CheckThreefold and CheckRemoved to FALSE when OptimizeBetween is TRUE
  if (OptimizeBetween) {
    CheckRemoved <- FALSE
    CheckThreefold <- FALSE
  }

  # Create empty data frame
  dfx <- data.frame(Data)
  dfxempty <- dfx[0,]

  # Check if correct input is provided and library can be run
  if (!is.numeric(as.matrix(na.omit(dfx)))) {
    stop("opdisDownsampling: Only numeric data allowed. Nothing to downsample.")
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
    warning(sprintf("opdisDownsampling: Size (%d) >= number of rows (%d). Nothing to downsample.",
                    Size, nrow(dfx)), call. = FALSE)
    return(list(ReducedData = dfx, RemovedData = dfxempty, ReducedInstances = rownames(dfx)))
  }

  if (Size <= 0) {
    warning(sprintf("opdisDownsampling: Size (%d) <= 0. All data will be removed.", Size), call. = FALSE)
    return(list(ReducedData = dfxempty, RemovedData = dfx, ReducedInstances = rownames(dfxempty)))
  }

  # Handle test statistic
  TestStats <- c("ad", "kuiper", "cvm", "wass", "dts", "ks", "kld", "amrdd", "euc", "nent")
  if (!(TestStat %in% TestStats)) {
    warning(paste0("opdisDownsampling: Possible TestStat = ", paste(TestStats, collapse = ", "),
                   ". TestStat set to default = 'ad'."), call. = FALSE)
    TestStat <- "ad"
  }

  # Initialize environment
  if (missing(Seed)) {
    Seed <- as.integer(get_seed()[1])
  }
  list.of.seeds <- 1:nTrials + Seed - 1
  nProc <- determine_n_cores(MaxCores)

  # Determine optimal chunk size if not provided or handle no-chunking request
  if (JobSize <= 0) {
    # 0 or negative means no chunking
    JobSize <- nTrials # Process everything in one chunk

    # Still print diagnostics if requested, but indicate no chunking
    if (verbose) {
      data_size_mb <- (nrow(dfx) * (ncol(dfx) - 1) * 8) / (1024 ^ 2)
      message(sprintf("Chunk size diagnostics:"))
      message(sprintf("  Data: %d rows x %d cols (%.1f MB)", nrow(dfx), ncol(dfx) - 1, data_size_mb))
      message(sprintf("  Trials: %d, Chunk size: %d (no chunking - single batch)",
                      nTrials, JobSize))
    }
  } else {
    # Use the provided JobSize or calculate optimal size
    if (is.null(JobSize)) {
      # This condition won't be hit with the new default, but kept for clarity
      JobSize <- calculate_optimal_chunk_size(
        n_rows = nrow(dfx),
        n_cols = ncol(dfx) - 1,
        nTrials = nTrials,
        nProc = nProc
      )
    }

    # Print diagnostics if requested
    print_chunk_diagnostics(nrow(dfx), ncol(dfx) - 1, nTrials, JobSize, verbose)
  }

  # Perform sampling and analyze picked data subsets
  ReducedDiag <- sample_and_analyze(DataAndClasses = dfx,
                                    TestStat = TestStat,
                                    Size = Size,
                                    list.of.seeds = list.of.seeds,
                                    PCAimportance = PCAimportance,
                                    NonNoiseSelection = NonNoiseSelection,
                                    UniformTestStat = UniformTestStat,
                                    UniformThreshold = UniformThreshold,
                                    nProc = nProc,
                                    CheckRemoved = CheckRemoved,
                                    CheckThreefold = CheckThreefold,
                                    OptimizeBetween = OptimizeBetween,
                                    JobSize = JobSize)

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
  tryCatch({
    AD_reduced_statMat <- matrix(NA_real_, nrow = n_trials, ncol = n_vars,
                                 dimnames = list(NULL, var_names))
    AD_removed_statMat <- matrix(NA_real_, nrow = n_trials, ncol = n_vars,
                                 dimnames = list(NULL, var_names))
    AD_reduced_vs_removed_statMat <- matrix(NA_real_, nrow = n_trials, ncol = n_vars,
                                            dimnames = list(NULL, var_names))
  }, error = function(e) {
    stop(sprintf("opdisDownsampling: Failed to allocate matrices: %s", e$message))
  })

  # Fill matrices with enhanced validation
  for (i in seq_len(n_trials)) {
    if (!is.null(ReducedDiag[[i]])) {
      tryCatch({
        # Validate data structure for each trial
        if (is.list(ReducedDiag[[i]]) &&
          all(c("ADv_reduced", "ADv_removed", "ADv_reduced_vs_removed") %in% names(ReducedDiag[[i]]))) {

          AD_reduced_statMat[i,] <- ReducedDiag[[i]][["ADv_reduced"]]
          AD_removed_statMat[i,] <- ReducedDiag[[i]][["ADv_removed"]]
          AD_reduced_vs_removed_statMat[i,] <- ReducedDiag[[i]][["ADv_reduced_vs_removed"]]
        } else {
          warning(sprintf("opdisDownsampling: Invalid data structure in trial %d, skipping.", i),
                  call. = FALSE)
        }
      }, error = function(e) {
        warning(sprintf("opdisDownsampling: Error processing trial %d: %s", i, e$message),
                call. = FALSE)
      })
    }
  }

  # Clean up the large ReducedDiag object early
  rm(ReducedDiag)
  gc()

  # Validate matrices before proceeding
  validate_matrices(list(AD_reduced_statMat, AD_removed_statMat, AD_reduced_vs_removed_statMat),
                    c("AD_reduced_statMat", "AD_removed_statMat", "AD_reduced_vs_removed_statMat"))

  # Find best subsample using refactored selection logic with error handling
  BestTrial <- tryCatch({
    if (OptimizeBetween) {
      select_best_trial_optimize_between(AD_reduced_vs_removed_statMat, WorstSample)
    } else if (CheckThreefold && CheckRemoved) {
      select_best_trial_threefold(AD_reduced_statMat, AD_removed_statMat, AD_reduced_vs_removed_statMat,WorstSample)
    } else if (CheckRemoved) {
      select_best_trial_check_removed(AD_reduced_statMat, AD_removed_statMat, WorstSample)
    } else {
      select_best_trial_reduced_only(AD_reduced_statMat, WorstSample)
    }
  }, error = function(e) {
    warning(sprintf("opdisDownsampling: Error in trial selection: %s. Using first trial.", e$message),
            call. = FALSE)
    1
  })

  # Validate BestTrial result
  if (is.na(BestTrial) || BestTrial < 1 || BestTrial > n_trials) {
    warning(sprintf("opdisDownsampling: Invalid best trial index (%s). Using first trial.",
                    as.character(BestTrial)), call. = FALSE)
    BestTrial <- 1
  }

  # Clean up stat matrices - moved after BestTrial validation
  rm(AD_reduced_statMat, AD_removed_statMat, AD_reduced_vs_removed_statMat)
  gc()

  # Get the best subsample
  df_reduced_final <- MakeReducedDataMat(DataAndClasses = dfx, Size = Size, Seed = list.of.seeds[BestTrial])
  ReducedData <- df_reduced_final$ReducedDataList
  RemovedData <- df_reduced_final$RemovedDataList

  if (Clsarg_missing) {
    ReducedData <- ReducedData[1:(ncol(ReducedData) - 1)]
    RemovedData <- RemovedData[1:(ncol(RemovedData) - 1)]
  }

  return(list(ReducedData = ReducedData,
              RemovedData = RemovedData,
              ReducedInstances = rownames(ReducedData),
              RemovedInstances = rownames(RemovedData)))
}