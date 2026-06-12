#' Perform Random Class-Proportional Downsampling
#'
#' This internal helper performs random class-proportional downsampling for a
#' given seed and returns both the selected and unselected rows.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param Size Desired size of the downsampled dataset, after conversion to an
#'   absolute number of rows by \code{opdisDownsampling()}.
#' @param Seed An integer value used as the random seed.
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item \code{ReducedDataList}: The downsampled dataset.
#'     \item \code{RemovedDataList}: The data not selected for the downsampled dataset.
#'   }
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

  # Validate that no variable in the reduced data contains only NAs
  # (excluding the Cls column)
  data_cols <- setdiff(names(ReducedDataList), "Cls")

  for (col in data_cols) {
    if (all(is.na(ReducedDataList[[col]]))) {
      warning(sprintf(
        "opdisDownsampling: Variable '%s' in reduced data contains only NA values with seed %d. This subsample may not be suitable.",
        col, Seed
      ), call. = FALSE)
    }
  }

  # Also check removed data if it will be used
  for (col in data_cols) {
    if (all(is.na(RemovedDataList[[col]]))) {
      warning(sprintf(
        "opdisDownsampling: Variable '%s' in removed data contains only NA values with seed %d. This subsample may not be suitable.",
        col, Seed
      ), call. = FALSE)
    }
  }

  return(list(
    ReducedDataList = ReducedDataList,
    RemovedDataList = RemovedDataList
  ))
}
