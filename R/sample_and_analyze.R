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
#'
#' @return A list of the results from the `make_and_analyse_subsample` function for
#'   each seed in `list.of.seeds`.
#'
#' @importFrom stats prcomp
#' @importFrom pbmcapply pbmclapply
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
sample_and_analyze <- function(DataAndClasses, TestStat, Size, list.of.seeds, PCAimportance, nProc, CheckRemoved) {
  # Identify relevant variables according to PCA projection, if selected
  selectedVars <- names(DataAndClasses)[1:(ncol(DataAndClasses)-1)]
  if (PCAimportance && length(list.of.seeds) > 1 && ncol(DataAndClasses) > 2) {
    pca1 <- prcomp(DataAndClasses[1:(ncol(DataAndClasses)-1)], retx = TRUE, center = TRUE, scale = TRUE)
    selectedVars_pca <- names(DataAndClasses)[which(names(DataAndClasses) %in% relevant_PCAvariables(pca1))]
    if (length(selectedVars_pca) > 0) selectedVars <- selectedVars_pca
  }
  
  # Central sampling function
  make_and_analyse_subsample <- function(DataAndClasses, TestStat, Size, Seed) {
    df_reduced <- MakeReducedDataMat(DataAndClasses, Size, Seed)
    ADv_reduced <- CompareReducedDataMat(DataAndClasses = DataAndClasses, ReducedDataList = df_reduced$ReducedDataList, TestStat = TestStat)
    if (CheckRemoved) {
      ADv_removed <- CompareReducedDataMat(DataAndClasses = DataAndClasses, ReducedDataList = df_reduced$RemovedDataList, TestStat = TestStat)  
    } else 
    {
      ADv_removed <- ADv_reduced
      ADv_removed[selectedVars] <- NA
    }
    return(list(ADv_reduced = ADv_reduced[selectedVars], ADv_removed = ADv_removed[selectedVars]) )
  }
  
  # Perform sampling and analyze picked data subsets
  if (nProc > 1) {
    if (Sys.info()[["sysname"]] == "Windows") {
      requireNamespace("foreach")
      doParallel::registerDoParallel(nProc)
      seed <- integer()
      ReducedDataMat <- foreach::foreach(seed = list.of.seeds) %dopar% {
        make_and_analyse_subsample(DataAndClasses, TestStat, Size, seed)
      }
      doParallel::stopImplicitCluster()
    } else {
      ReducedDataMat <- pbmcapply::pbmclapply(
        list.of.seeds,
        function(seed) make_and_analyse_subsample(DataAndClasses, TestStat, Size, seed),
        mc.cores = nProc
      )
    }
  } else {
    ReducedDataMat <- lapply(
      list.of.seeds,
      function(seed) make_and_analyse_subsample(DataAndClasses, TestStat, Size, seed)
    )
  }
  
  return(ReducedDataMat)
}
