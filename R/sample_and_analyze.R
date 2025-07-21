#' Perform Sampling and Analyze Picked Data Subsets
#'
#' This function performs sampling and analyzes the picked data subsets using the
#' `make_and_analyse_subsample` function.
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
#' @param nProc The number of cores to use for parallel processing.
#' @param CheckRemoved A logical value indicating whether to also optimize the removed part
#'   of the data for distribution equality with the original.
#' @param JobSize Number of seeds to process in each chunk (default: min(50, length(list.of.seeds)/nProc))
#'
#' @return A list of the results from the `make_and_analyse_subsample` function for
#'   each seed in `list.of.seeds`.
#'
#' @importFrom stats prcomp
#' @importFrom pbmcapply pbmclapply
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
sample_and_analyze <- function(DataAndClasses, TestStat, Size, list.of.seeds, PCAimportance, nProc, CheckRemoved, JobSize = NULL) {

  # Set default chunk size if not provided
  if (is.null(JobSize)) {
    JobSize <- min(50, max(1, ceiling(length(list.of.seeds) / max(1, nProc))))
  }

  # Identify relevant variables according to PCA projection, if selected
  selectedVars <- names(DataAndClasses)[1:(ncol(DataAndClasses)-1)]
  if (PCAimportance && length(list.of.seeds) > 1 && ncol(DataAndClasses) > 2) {
    # Use less memory-intensive PCA options
    pca1 <- prcomp(DataAndClasses[1:(ncol(DataAndClasses)-1)],
                   retx = FALSE,  # Don't store transformed data
                   center = TRUE,
                   scale = TRUE,
                   rank. = min(10, ncol(DataAndClasses)-1))  # Limit components
    selectedVars_pca <- names(DataAndClasses)[which(names(DataAndClasses) %in% relevant_PCAvariables(pca1))]
    if (length(selectedVars_pca) > 0) selectedVars <- selectedVars_pca
    rm(pca1)  # Clean up PCA object
    gc()      # Force garbage collection
  }

  # Pre-compute data subset for selected variables only
  DataSubset <- DataAndClasses[, c(selectedVars, names(DataAndClasses)[ncol(DataAndClasses)]), drop = FALSE]

  # Memory-optimized sampling function
  make_and_analyse_subsample <- function(DataSubset, TestStat, Size, Seed, selectedVars, CheckRemoved) {
    df_reduced <- MakeReducedDataMat(DataSubset, Size, Seed)
    ADv_reduced <- CompareReducedDataMat(DataAndClasses = DataSubset,
                                         ReducedDataList = df_reduced$ReducedDataList,
                                         TestStat = TestStat)

    if (CheckRemoved) {
      ADv_removed <- CompareReducedDataMat(DataAndClasses = DataSubset,
                                           ReducedDataList = df_reduced$RemovedDataList,
                                           TestStat = TestStat)
    } else {
      ADv_removed <- rep(NA_real_, length(selectedVars))
      names(ADv_removed) <- selectedVars
    }

    # Clean up intermediate objects
    rm(df_reduced)

    return(list(ADv_reduced = ADv_reduced[selectedVars],
                ADv_removed = ADv_removed[selectedVars]))
  }

  # Process in chunks to reduce memory pressure
  process_chunk <- function(seed_chunk, DataSubset, TestStat, Size, selectedVars, CheckRemoved, use_parallel = FALSE, cores = 1) {
    if (use_parallel && cores > 1 && length(seed_chunk) > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        # Use foreach for Windows
        doParallel::registerDoParallel(min(cores, length(seed_chunk)))
        seed <- integer()
        result <- foreach::foreach(seed = seed_chunk,
                                   .combine = 'c',
                                   .maxcombine = length(seed_chunk),
                                   .multicombine = TRUE,
                                   .packages = c()) %dopar% {
          list(make_and_analyse_subsample(DataSubset, TestStat, Size, seed, selectedVars, CheckRemoved))
        }
        doParallel::stopImplicitCluster()
      } else {
        # Use mclapply for Unix-like systems
        result <- pbmcapply::pbmclapply(
          seed_chunk,
          function(seed) make_and_analyse_subsample(DataSubset, TestStat, Size, seed, selectedVars, CheckRemoved),
          mc.cores = min(cores, length(seed_chunk)),
          mc.preschedule = TRUE  # Better load balancing
        )
      }
    } else {
      # Sequential processing
      result <- lapply_with_bar(
        seed_chunk,
        function(seed) make_and_analyse_subsample(DataSubset, TestStat, Size, seed, selectedVars, CheckRemoved)
      )
    }
    return(result)
  }

  # Split seeds into chunks
  seed_chunks <- split(list.of.seeds, ceiling(seq_along(list.of.seeds) / JobSize))

  # Initialize result list
  ReducedDataMat <- vector("list", length(list.of.seeds))
  current_index <- 1

  # Process chunks
  for (chunk in seed_chunks) {
    chunk_result <- process_chunk(
      seed_chunk = chunk,
      DataSubset = DataSubset,
      TestStat = TestStat,
      Size = Size,
      selectedVars = selectedVars,
      CheckRemoved = CheckRemoved,
      use_parallel = nProc > 1,
      cores = nProc
    )

    # Store results
    end_index <- current_index + length(chunk) - 1
    ReducedDataMat[current_index:end_index] <- chunk_result
    current_index <- end_index + 1

    # Force garbage collection after each chunk
    gc()
  }

  return(ReducedDataMat)
}
