# ===============================================================================
# SUBSAMPLE ANALYSIS FUNCTIONS
# ===============================================================================
#
# This file contains the core function for creating and analyzing individual
# subsamples during the optimization process.

#' Create and Analyze a Single Subsample
#'
#' This function creates a single reduced dataset using a specific seed and
#' analyzes its distributional similarity to the original data and/or removed data.
#' This is the core function called for each trial during the optimization process.
#'
#' @param DataSubset A data frame containing the data subset to be analyzed
#' @param TestStat Statistical test for comparing distributions
#' @param Size Desired size of the reduced dataset
#' @param Seed Random seed for this specific trial
#' @param selectedVars Character vector of variable names to include
#' @param CheckRemoved Whether to compute removed vs original comparison
#' @param CheckThreefold Whether to compute reduced vs removed comparison
#' @param OptimizeBetween Whether to focus only on reduced vs removed optimization
#' @return Named list with ADv_reduced, ADv_removed, ADv_reduced_vs_removed vectors
#' @keywords internal
make_and_analyse_subsample <- function(DataSubset, TestStat, Size, Seed, selectedVars, CheckRemoved, CheckThreefold, OptimizeBetween) {

  # Create the data split
  df_reduced <- MakeReducedDataMat(DataSubset, Size, Seed)

  # Always compute reduced vs original comparison
  ADv_reduced <- CompareReducedDataMat(DataAndClasses = DataSubset,
                                       ReducedDataList = df_reduced$ReducedDataList,
                                       TestStat = TestStat)

  # Compute removed vs original comparison if needed
  if (CheckRemoved || OptimizeBetween) {
    ADv_removed <- CompareReducedDataMat(DataAndClasses = DataSubset,
                                         ReducedDataList = df_reduced$RemovedDataList,
                                         TestStat = TestStat)
  } else {
    ADv_removed <- rep(NA_real_, length(selectedVars))
    names(ADv_removed) <- selectedVars
  }

  # Compute reduced vs removed comparison if needed
  if (CheckThreefold || OptimizeBetween) {
    ADv_reduced_vs_removed <- CompareReducedDataMat(DataAndClasses = df_reduced$ReducedDataList,
                                                    ReducedDataList = df_reduced$RemovedDataList,
                                                    TestStat = TestStat)
  } else {
    ADv_reduced_vs_removed <- rep(NA_real_, length(selectedVars))
    names(ADv_reduced_vs_removed) <- selectedVars
  }

  # Clean up intermediate objects
  rm(df_reduced)

  return(list(ADv_reduced = ADv_reduced[selectedVars],
              ADv_removed = ADv_removed[selectedVars],
              ADv_reduced_vs_removed = ADv_reduced_vs_removed[selectedVars]))
}