#' Compare the Distributions of Original and Downsampled Datasets
#'
#' This function compares the distributions of the original and downsampled datasets
#' using a specified statistical test.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param ReducedDataList The downsampled dataset.
#' @param TestStat A character string specifying the statistical test to be used for
#'   comparing the distributions.
#'
#' @return The values of the specified statistical test comparing the original and
#'   downsampled distributions.
#'
CompareReducedDataMat <- function(DataAndClasses, ReducedDataList, TestStat) {
  # Compare the original and downsampled distributions
  return(mapply(
    CompDistrib,
    vector1 = DataAndClasses[ 1:(ncol(DataAndClasses) - 1)],
    vector2 = ReducedDataList[ 1:(ncol(DataAndClasses) - 1)],
    MoreArgs = list(TestStat = TestStat)
  ))
}
