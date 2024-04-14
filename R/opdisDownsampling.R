#' Perform Class-Proportional Downsampling
#'
#' This function samples a subset of data based on the similarity of its probability
#' distribution to that of the original data.
#'
#' @param Data A data frame containing the data to be downsampled.
#' @param Cls A vector of class labels for the data. If not provided, all instances
#'   will be assigned to a single class.
#' @param Size A numeric value between 0 and 1 representing the desired size of the
#'   downsampled dataset as a proportion of the original dataset.
#' @param Seed An integer value to be used as the random seed for reproducibility.
#' @param nTrials The number of trials to perform when sampling the data.
#' @param TestStat A character string specifying the statistical test to be used for
#'   comparing the distributions.
#' @param MaxCores The maximum number of cores to use for parallel processing.
#' @param PCAimportance A logical value indicating whether to use PCA to identify
#'   relevant variables.
#'
#' @return A list with the following elements:
#'   - `ReducedData`: The downsampled dataset.
#'   - `RemovedData`: The removed data from the original dataset.
#'   - `ReducedInstances`: The row names of the downsampled dataset.
#'
#' @useDynLib opdisDownsampling, .registration = TRUE
#' @importFrom methods hasArg
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom pbmcapply pbmclapply
#' @import foreach
#' @export
opdisDownsampling <- function(Data, Cls, Size, Seed, nTrials = 1000, TestStat = "ad",
                              MaxCores = getOption("mc.cores", 2L), PCAimportance = FALSE) {
  dfx <- data.frame(Data)
  dfxempty <- dfx[0, ]
  
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
  
  # Handle size parameter
  if (Size >= nrow(dfx)) {
    warning("opdisDownsampling: Size >= length of 'Data'. Nothing to downsample.", call. = FALSE)
    return(list(ReducedData = dfx, RemovedData = dfxempty, ReducedInstances = rownames(dfx)))
  }
  
  if (Size <= 0) {
    warning("opdisDownsampling: Size <= 0. All data will be removed.", call. = FALSE)
    return(list(ReducedData = dfxempty, RemovedData = dfx, ReducedInstances = rownames(dfxempty)))
  }
  
  # Handle test statistic
  TestStats <- c("ad", "kuiper", "cvm", "wass", "dts", "ks", "kld", "amrdd", "euc")
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
  num_workers <- parallel::detectCores()
  nProc <- min(num_workers - 1, MaxCores)
  
  # Perform sampling and analyze picked data subsets
  ReducedDiag <- sample_and_analyze(DataAndClasses = dfx,
                                    TestStat = TestStat,
                                    Size = Size,
                                    list.of.seeds = list.of.seeds,
                                    PCAimportance = PCAimportance,
                                    nProc = nProc)
  
  # Find best subsample
  ADstatMat <- do.call(rbind, ReducedDiag)
  BestTrial <- which.min(apply(ADstatMat, 1, max))
  
  # Get the best subsample
  df_reduced_final <- MakeReducedDataMat(DataAndClasses = dfx, Size = Size, Seed = list.of.seeds[BestTrial])
  ReducedData <- df_reduced_final$ReducedDataList
  RemovedData <- df_reduced_final$RemovedDataList
  
  if (Clsarg_missing) {
    ReducedData <- ReducedData[1:(ncol(ReducedData) - 1)]
    RemovedData <- RemovedData[1:(ncol(RemovedData) - 1)]
  }
  
  return(list(ReducedData = ReducedData, RemovedData = RemovedData, ReducedInstances = rownames(ReducedData)))
}
