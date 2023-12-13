# Performs random class-proportional downsampling and compares the sampled and
# original distributions
#' @importFrom caTools sample.split
MakeReducedDataMat <-
  function(DataAndClasses, TestStat, Size, Seed) {
    set.seed(Seed)
    sample <-
      caTools::sample.split(DataAndClasses$Cls, SplitRatio = Size / nrow(DataAndClasses))
    ReducedDataList <- subset(DataAndClasses, sample == TRUE)
    RemovedDataList <- subset(DataAndClasses, sample == FALSE)
    ADv <-
      mapply(CompDistrib,
             vector1 = DataAndClasses[1:(ncol(DataAndClasses) - 1)],
             vector2 = ReducedDataList[1:(ncol(ReducedDataList) -
                                            1)],
             TestStat)
    return(
      list(
        ReducedDataList = ReducedDataList,
        RemovedDataList = RemovedDataList,
        ADv = ADv
      )
    )
  }