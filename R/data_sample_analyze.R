#' Perform Sampling and Analyze Picked Data Subsets
#'
#' This function performs sampling and analyzes the picked data subsets using the
#' `make_and_analyse_subsample` function with modular variable selection and
#' chunk-based processing for memory efficiency.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param TestStat A character string specifying the statistical test to be used for
#'   comparing the distributions.
#' @param Size A numeric value between 0 and 1 representing the desired size of the
#'   downsampled dataset as a proportion of the original dataset.
#' @param list.of.seeds A vector of integer values to be used as the random seeds for
#'   reproducibility.
#' @param PCAimportance A logical value indicating whether to use PCA to identify
#'   relevant variables.
#' @param NonNoiseSelection A logical value indicating whether to use non-uniform
#'   distribution tests to identify relevant variables.
#' @param UniformTestStat A character string specifying the statistical test to be used for
#'   non-uniform variable selection. Available options are: "ks" (Kolmogorov-Smirnov),
#'   "ad" (Anderson-Darling), "kuiper", "cvm" (Cramér-von Mises), "wass" (Wasserstein),
#'   "dts" (Distributional Transform Statistic), "kld" (Kullback-Leibler divergence),
#'   "amrdd" (Average Mean Root of Distributional Differences), and "euc" (Euclidean distance).
#'   Only used when NonNoiseSelection = TRUE (default: "ks").
#' @param UniformThreshold Threshold value for non-uniform variable selection (default: 0.1).
#' @param nProc The number of cores to use for parallel processing.
#' @param CheckRemoved A logical value indicating whether to also optimize the removed part
#'   of the data for distribution equality with the original.
#' @param CheckThreefold A logical value indicating whether to also optimize the reduced part
#'   of the data for distribution equality with the removed part. Ignored when CheckRemoved is FALSE.
#' @param OptimizeBetween A logical value indicating whether to optimize the reduced part
#'   of the data for distribution equality with the removed part. If set, all other comparisons are not performed.
#' @param JobSize Number of seeds to process in each chunk (default: min(50, length(list.of.seeds)/nProc))
#'
#' @return A list of the results from the `make_and_analyse_subsample` function for
#'   each seed in `list.of.seeds`.
#'
#' @details The function operates in the following phases:
#' \itemize{
#'   \item Variable Selection: Applies PCA and/or non-uniform distribution tests to identify relevant variables
#'   \item Data Preparation: Creates a subset containing only selected variables and class labels
#'   \item Chunk Processing: Processes seeds in chunks to manage memory efficiently
#'   \item Parallel Execution: Uses platform-appropriate parallel processing when beneficial
#' }
#'
#' @importFrom stats prcomp
#' @importFrom pbmcapply pbmclapply
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
sample_and_analyze <- function(DataAndClasses, TestStat, Size, list.of.seeds, PCAimportance, nProc,
                               CheckRemoved, CheckThreefold, OptimizeBetween, JobSize = NULL,
                               NonNoiseSelection = FALSE, UniformTestStat = "ks",
                               UniformThreshold = 0.05) {

  # Set default chunk size if not provided
  if (is.null(JobSize)) {
    JobSize <- min(50, max(1, ceiling(length(list.of.seeds) / max(1, nProc))))
  }

  # Phase 1: Variable Selection
  # ==========================
  selectedVars <- select_variables(
    DataAndClasses = DataAndClasses,
    PCAimportance = PCAimportance,
    NonNoiseSelection = NonNoiseSelection,
    UniformTestStat = UniformTestStat,
    UniformThreshold = UniformThreshold,
    list.of.seeds = list.of.seeds
  )

  # Print variable selection summary for user feedback
  all_vars <- names(DataAndClasses)[1:(ncol(DataAndClasses) - 1)]

  if (PCAimportance || NonNoiseSelection) {
    cat(paste("Final variable selection:", length(selectedVars), "out of",
              length(all_vars), "variables\n"))
    cat(paste("Selected variables:", paste(selectedVars, collapse = ", "), "\n"))
  }

  # Phase 2: Data Preparation
  # =========================
  # Pre-compute data subset for selected variables only to reduce memory usage
  DataSubset <- DataAndClasses[, c(selectedVars, names(DataAndClasses)[ncol(DataAndClasses)]), drop = FALSE]

  # Prepare processing parameters for chunk processing
  processing_params <- list(
    TestStat = TestStat,
    Size = Size,
    selectedVars = selectedVars,
    CheckRemoved = CheckRemoved,
    CheckThreefold = CheckThreefold,
    OptimizeBetween = OptimizeBetween
  )

  # Phase 3: Chunk-Based Processing
  # ===============================
  # Process seeds in memory-efficient chunks with appropriate parallelization
  ReducedDataMat <- process_seeds_in_chunks(
    list.of.seeds = list.of.seeds,
    DataSubset = DataSubset,
    processing_params = processing_params,
    JobSize = JobSize,
    nProc = nProc
  )

  return(ReducedDataMat)
}