#' Perform Sampling and Analyze Picked Data Subsets
#'
#' This internal helper selects the variables used for distribution comparison,
#' prepares the reduced working data set, and processes candidate sampling seeds
#' in chunks.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param TestStat A character string specifying the statistical test to be used for
#'   comparing the distributions.
#' @param Size Desired size of the downsampled dataset, after conversion to an
#'   absolute number of rows by \code{opdisDownsampling()}.
#' @param list.of.seeds A vector of integer values to be used as random seeds.
#' @param PCAimportance A logical value indicating whether to use PCA-based
#'   variable selection.
#' @param nProc The number of cores to use for parallel processing.
#' @param JobSize Number of seeds to process in each chunk.
#'
#' @return A list of results from \code{make_and_analyse_subsample()}, one result
#'   for each seed in \code{list.of.seeds}.
#'
#' @details The function operates in the following phases:
#' \itemize{
#'   \item optional PCA-based variable selection,
#'   \item data preparation using the selected variables and class labels,
#'   \item chunked seed processing to control memory use,
#'   \item parallel or sequential execution depending on \code{nProc}.
#' }
#'
#' @importFrom stats prcomp
#' @importFrom pbmcapply pbmclapply
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
sample_and_analyze <- function(DataAndClasses, TestStat, Size, list.of.seeds, PCAimportance, nProc,
                               JobSize = NULL) {
  # Set default chunk size if not provided
  if (is.null(JobSize)) {
    JobSize <- min(50, max(1, ceiling(length(list.of.seeds) / max(1, nProc))))
  }

  # Phase 1: Variable Selection
  selectedVars <- select_variables(
    DataAndClasses = DataAndClasses,
    PCAimportance = PCAimportance,
    list.of.seeds = list.of.seeds
  )

  # Print variable selection summary for user feedback
  all_vars <- names(DataAndClasses)[1:(ncol(DataAndClasses) - 1)]

  if (PCAimportance) {
    cat(paste(
      "Final variable selection:", length(selectedVars), "out of",
      length(all_vars), "variables\n"
    ))
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
    selectedVars = selectedVars
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
