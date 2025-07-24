#' Perform Random Class-Proportional Downsampling and Compare Distributions
#'
#' This function performs random class-proportional downsampling on a dataset and
#' compares the sampled and original distributions using a specified statistical test.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param Size A numeric value between 0 and 1 representing the desired size of the
#'   downsampled dataset as a proportion of the original dataset.
#' @param Seed An integer value to be used as the random seed for reproducibility.
#'
#' @return A list with the following elements:
#'   - `ReducedDataList`: The downsampled dataset.
#'   - `RemovedDataList`: The removed data from the original dataset.
#'
#' @importFrom caTools sample.split
#'
MakeReducedDataMat <- function(DataAndClasses, Size, Seed) {
  # Set the random seed for reproducibility
  set.seed(Seed)
  
  # Perform class-proportional downsampling
  sample <- caTools::sample.split(DataAndClasses$Cls, SplitRatio = Size)
  ReducedDataList <- DataAndClasses[sample, ]
  RemovedDataList <- DataAndClasses[!sample, ]
  
  return(list(
    ReducedDataList = ReducedDataList,
    RemovedDataList = RemovedDataList
  ))
}


