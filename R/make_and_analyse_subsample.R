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
#' @return Named list with ADv_reduced, ADv_removed
#' @keywords internal
make_and_analyse_subsample <- function(DataSubset, TestStat, Size, Seed, selectedVars) {
  # Create the data split
  df_reduced <- MakeReducedDataMat(DataSubset, Size, Seed)

  # Always compute reduced vs original comparison
  ADv_reduced <- CompareReducedDataMat(
    DataAndClasses = DataSubset,
    ReducedDataList = df_reduced$ReducedDataList,
    TestStat = TestStat
  )

  # Compute removed vs original comparison if needed
  ADv_removed <- rep(NA_real_, length(selectedVars))
  names(ADv_removed) <- selectedVars

  # Clean up intermediate objects
  rm(df_reduced)

  return(list(
    ADv_reduced = ADv_reduced[selectedVars],
    ADv_removed = ADv_removed[selectedVars]
  ))
}
