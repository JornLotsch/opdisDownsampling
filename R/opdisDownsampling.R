#Calculates the Impact effect size measure.
#Samples a subset of data based on the similarity of
#its probability distribution to that of the original data
#' @useDynLib(opdisDownsampling, .registration = TRUE)
#' @importFrom caTools sample.split
#' @importFrom methods hasArg
#' @importFrom twosamples ad_stat kuiper_stat cvm_stat wass_stat dts_stat
#' @importFrom stats density dist ks.test prcomp
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @importFrom EucDist EucDist
#' @importFrom KullbLeiblKLD2 KullbLeiblKLD2
#' @export
opdisDownsampling <- function(Data, Cls, Size, Seed, nTrials, TestStat = "ad",
                             MaxCores = 8, JobSize = 10000) {
  dfx <- data.frame(Data)
  if (hasArg("Cls") == TRUE) {
    if (length(Cls) != nrow(dfx)) {
      stop("opdisDownsampling: Unequal number of cases and class memberships.")
    }
    dfx$Cls <- Cls
  } else {
    dfx$Cls <- 1
  }
  if (Size >= nrow(dfx)) {
    warning("opdisDownsampling: Size >= length of 'Data'. Nothing to downsample.", call. = FALSE)
    ReducedData <- dfx
    RemovedData <- vector()
  } else {
    if (hasArg("nTrials") == TRUE) nTrials = nTrials else nTrials = 1000
    if (hasArg("Seed") == TRUE) Seed = Seed else Seed = 42
    if (hasArg("TestStat") == FALSE) TestStat = "ad" else TestStat = TestStat

    requireNamespace("twosamples")
    CompDistrib <- function(vector1, vector2) {
      switch(TestStat,
             "ad"  = {ad_stat(vector1, vector2)},
             "kuiper" = {kuiper_stat(vector1, vector2)},
             "cvm" = {cvm_stat(vector1, vector2)},
             "wass" = {wass_stat(vector1, vector2)},
             "dts" = {dts_stat(vector1, vector2)},
             "ks" = {ks.test(vector1, vector2)$statistic},
             "kld" = {KullbLeiblKLD2(vector1, vector2)$KLD},
             "amrdd" = {amrdd(vector1, vector2)},
             "euc" = {EucDist(vector1, vector2)}
      )
    }

    list.of.seeds.all <- 1:nTrials  + Seed
    list.of.seeds <- split(list.of.seeds.all, ceiling(seq_along(list.of.seeds.all) / JobSize))
    nlist.of.seeds <- unlist(lapply(list.of.seeds, length))
    ADstatAll <- vector()
    ReducedDataI <- list()
    RemovedDataI <- list()
    for (i in 1:length(list.of.seeds)) {
      ADstat <- vector()
      requireNamespace("caTools")
      if (nlist.of.seeds[i] * length(Cls) > 6000 & .Platform$OS.type != "windows" & MaxCores > 1) {
        requireNamespace("parallel")
        chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
        if (nzchar(chk) && chk == "TRUE") {
          # use 2 cores in CRAN/Travis/AppVeyor
          num_workers <- 2L
        } else {
          # use all cores in devtools::test()
          num_workers <- detectCores()
        }
        nProc <- min(num_workers - 1, MaxCores)
        requireNamespace("pbmcapply")
        ReducedDataMat <-
          pbmclapply(1:nlist.of.seeds[i], function(x) {
            set.seed(list.of.seeds[[i]][x])
            sample <- sample.split(dfx$Cls, SplitRatio = Size / nrow(dfx))
            ReducedDataList <- subset(dfx, sample == TRUE)
            RemovedDataList <- subset(dfx, sample == FALSE)
            ADv <- mapply(CompDistrib, dfx[1:(ncol(dfx) - 1)], ReducedDataList[1:(ncol(ReducedDataList) - 1)])
            return(list(ReducedDataList = ReducedDataList, RemovedDataList = RemovedDataList, ADv = ADv))
          }, mc.cores = nProc)
      } else {
        ReducedDataMat <- lapply(1:nlist.of.seeds[i], function(x) {
          pb <- txtProgressBar(min = 0, max = nlist.of.seeds[i], style = 3)
          set.seed(list.of.seeds[[i]][x])
          sample <- sample.split(dfx$Cls, SplitRatio = Size / nrow(dfx))
          ReducedDataList <- subset(dfx, sample == TRUE)
          RemovedDataList <- subset(dfx, sample == FALSE)
          ADv <- mapply(CompDistrib, dfx[1:(ncol(dfx) - 1)], ReducedDataList[1:(ncol(ReducedDataList) - 1)])
          setTxtProgressBar(pb, x)
          return(list(ReducedDataList = ReducedDataList, RemovedDataList = RemovedDataList, ADv = ADv))
        })
      }
      ADstat <- rbind(ADstat, unlist(lapply(ReducedDataMat, "[[", "ADv")))
      ADstatMat <- data.frame(matrix(ADstat, ncol = nlist.of.seeds[i]))

      BestTrial <- which.min(apply(ADstatMat, 2, function(x) max(x)))
      BestTrialStat <- min(apply(ADstatMat, 2, function(x) max(x)))
      ADstatAll <- append(ADstatAll, BestTrialStat)
      ReducedDataI[[i]] <- ReducedDataMat[[BestTrial]][["ReducedDataList"]]
      RemovedDataI[[i]] <- ReducedDataMat[[BestTrial]][["RemovedDataList"]]
    }

    BestTrialAll <- which.min(ADstatAll)
    ReducedData <- ReducedDataI[[BestTrialAll]]
    RemovedData <- RemovedDataI[[BestTrialAll]]

    if (hasArg("Cls") == FALSE) {
      ReducedData <- as.vector(ReducedData$Data)
      RemovedData <- as.vector(RemovedData$Data)
    }
  }
  return(list(ReducedData = ReducedData, RemovedData = RemovedData))
}
