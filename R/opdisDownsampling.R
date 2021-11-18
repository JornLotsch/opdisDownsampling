# Samples a subset of data based on the similarity of its probability
# distribution to that of the original data
#' @useDynLib(opdisDownsampling, .registration = TRUE)
#' @importFrom methods hasArg
#' @importFrom parallel detectCores
#' @importFrom benchmarkme get_ram
#' @importFrom memuse Sys.meminfo
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @import foreach
#' @export
opdisDownsampling <- function(Data, Cls, Size, Seed, nTrials = 1000, TestStat = "ad",
  MaxCores = getOption("mc.cores", 2L), JobSize = 10000, PCAimportance = FALSE) {
  dfx <- data.frame(Data)
  if (hasArg("Cls") == TRUE) {
    if (length(Cls) != nrow(dfx)) {
      stop("opdisDownsampling: Unequal number of cases and class memberships.")
    }
  } else {
    Cls <- rep(1, nrow(dfx))
  }
  dfx$Cls <- Cls

  if (Size >= nrow(dfx)) {
    warning("opdisDownsampling: Size >= length of 'Data'.
    Nothing to downsample.",
      call. = FALSE)
    ReducedData <- dfx
    RemovedData <- vector()
  } else {
    if (is.numeric(as.matrix(na.omit(Data))) == FALSE) {
      warning("opdisDownsampling: Only numeric data allowed.
    Nothing to downsample.",
        call. = FALSE)
      ReducedData <- dfx
      RemovedData <- vector()
    } else {
      if (!missing(Seed))
        Seed <- Seed else Seed <- 42

      num_workers <- parallel::detectCores()
      nProc <- min(num_workers - 1, MaxCores)

      list.of.seeds.all <- 1:nTrials + Seed

      if (!missing(JobSize))
        JobSize <- JobSize else {
        switch(Sys.info()[["sysname"]], Windows = {
          JobSize <- as.numeric(memuse::Sys.meminfo()$freeram) * 0.8 * nTrials *
          (dim(subset(dfx, select = -c(Cls)))[1] * dim(subset(dfx, select = -c(Cls)))[2])
        }, JobSize <- as.numeric(benchmarkme::get_ram()) * 0.8 * nTrials *
          (dim(subset(dfx, select = -c(Cls)))[1] * dim(subset(dfx, select = -c(Cls)))[2]))
      }

      if (nProc > 1) {
        list.of.seeds <- split(list.of.seeds.all, ceiling(seq_along(list.of.seeds.all)/max(nTrials/nProc,
          JobSize)))
      } else {
        list.of.seeds <- split(list.of.seeds.all, 1)
      }

      # Main part. For reasons of computing speed, three separate
      # versions are available, i.e., for Windows, Linux, and for single-core
      # processing.
      nlist.of.seeds <- as.integer(unlist(lapply(list.of.seeds, length)))
      ADstatAll <- vector()
      ReducedDataI <- list()
      RemovedDataI <- list()
      for (i in 1:length(list.of.seeds)) {
        ADstat <- vector()
        if (nProc > 1) {
          switch(Sys.info()[["sysname"]], Windows = {
          requireNamespace("foreach")
          doParallel::registerDoParallel(nProc)
          x <- integer()
          ReducedDataMat <- foreach::foreach(x = 1:nlist.of.seeds[i]) %dopar%
            {
            MakeReducedDataMat(DataAndClasses = dfx, TestStat = TestStat,
              Size = Size, Seed = list.of.seeds[[i]][x])
            }
          doParallel::stopImplicitCluster()
          }, {
          ReducedDataMat <- parallel::mclapply(1:nlist.of.seeds[i], function(x) {
            MakeReducedDataMat(DataAndClasses = dfx, TestStat = TestStat,
            Size = Size, Seed = list.of.seeds[[i]][x])
          }, mc.cores = nProc)
          })
        } else {
          ReducedDataMat <- lapply(1:nlist.of.seeds[i], function(x) {
          MakeReducedDataMat(DataAndClasses = dfx, TestStat = TestStat,
            Size = Size, Seed = list.of.seeds[[i]][x])
          })
        }

        ADstat <- rbind(ADstat, unlist(lapply(ReducedDataMat, "[[", "ADv")))
        ADstatMat <- data.frame(matrix(ADstat, ncol = nlist.of.seeds[i]))

        if (PCAimportance == TRUE & nlist.of.seeds[i] > 1 & ncol(dfx) > 2) {
          pca1 <- prcomp(dfx[1:(ncol(dfx) - 1)], retx = TRUE, center = TRUE,
          scale = TRUE)
          is.integer0 <- function(x) {
          is.integer(x) && length(x) == 0L
          }
          selectedVars <- which(names(dfx) %in% relevant_PCAvariables(res.pca = pca1))
          if (is.integer0(selectedVars) == FALSE)
          ADstatMat <- ADstatMat[c(selectedVars), ]
        }

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
  }
  return(list(ReducedData = ReducedData, RemovedData = RemovedData, ReducedInstances = rownames(ReducedData)))
}
