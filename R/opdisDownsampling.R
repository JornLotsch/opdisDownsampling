# Samples a subset of data based on the similarity of its probability
# distribution to that of the original data
#' @useDynLib(opdisDownsampling, .registration = TRUE)
#' @importFrom caTools sample.split
#' @importFrom methods hasArg
#' @importFrom twosamples ad_stat kuiper_stat cvm_stat wass_stat dts_stat
#' @importFrom stats ks.test prcomp na.omit
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster parLapply
#' @importFrom pbmcapply pbmclapply
#' @importFrom EucDist EucDist
#' @importFrom KullbLeiblKLD2 KullbLeiblKLD2
#' @importFrom benchmarkme get_ram
#' @export
opdisDownsampling <- function(Data, Cls, Size, Seed, nTrials = 1000, TestStat = "ad",
  MaxCores = 2048, JobSize = 1000, PCAimportance = FALSE) {
  dfx <- data.frame(Data)
  if (hasArg("Cls") == TRUE) {
    if (length(Cls) != nrow(dfx)) {
      stop("opdisDownsampling: Unequal number of cases and class memberships.")
    }
    dfx$Cls <- Cls
  } else {
    dfx$Cls <- 1
    Cls <- as.vector(dfx$Cls)
  }
  if (Size >= nrow(dfx)) {
    warning("opdisDownsampling: Size >= length of 'Data'.
    Nothing to downsample.",
      call. = FALSE)
    ReducedData <- dfx
    RemovedData <- vector()
  } else {
    if (!missing(Seed))
      Seed <- Seed else Seed <- 42
    if (!missing(TestStat))
      TestStat <- TestStat else TestStat <- "ad"

    CompDistrib <- function(vector1, vector2) {
      if (length(vector1[!is.na(vector1)]) * length(vector2[!is.na(vector2)]) ==
        0) {
        Stat <- 1e+27
      } else {
        Stat <- switch(TestStat, ad = {
          twosamples::ad_stat(na.omit(vector1), na.omit(vector2))
        }, kuiper = {
          twosamples::kuiper_stat(na.omit(vector1), na.omit(vector2))
        }, cvm = {
          twosamples::cvm_stat(na.omit(vector1), na.omit(vector2))
        }, wass = {
          twosamples::wass_stat(na.omit(vector1), na.omit(vector2))
        }, dts = {
          twosamples::dts_stat(na.omit(vector1), na.omit(vector2))
        }, ks = {
          ks.test(na.omit(vector1), na.omit(vector2))$statistic
        }, kld = {
          KullbLeiblKLD2(na.omit(vector1), na.omit(vector2))$KLD
        }, amrdd = {
          amrdd(na.omit(vector1), na.omit(vector2))
        }, euc = {
          EucDist(na.omit(vector1), na.omit(vector2))
        })
      }
      return(Stat)
    }

    requireNamespace("parallel")
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      num_workers <- 2L
    } else {
      num_workers <- parallel::detectCores()
    }
    nProc <- min(num_workers - 1, MaxCores)

    mclapply.hack <- function(...) {
      ## Create a cluster
      size.of.list <- length(list(...)[[1]])
      cl <- makeCluster(min(size.of.list, detectCores()))
      ## Find out the names of the loaded packages
      loaded.package.names <- c(sessionInfo()$basePkgs, names(sessionInfo()$otherPkgs))
      tryCatch({
        ## Copy over all of the objects within scope to all clusters.
        this.env <- environment()
        while (identical(this.env, globalenv()) == FALSE) {
          clusterExport(cl, ls(all.names = TRUE, envir = this.env), envir = this.env)
          this.env <- parent.env(environment())
        }
        clusterExport(cl, ls(all.names = TRUE, envir = globalenv()), envir = globalenv())

        ## Load the libraries on all the clusters N.B. length(cl) returns the
        ## number of clusters
        parLapply(cl, 1:length(cl), function(xx) {
          lapply(loaded.package.names, function(yy) {
          require(yy, character.only = TRUE)
          })
        })

        ## Run the lapply in parallel
        return(parLapply(cl, ...))
      }, finally = {
        ## Stop the cluster
        stopCluster(cl)
      })
    }

    mclapply <- switch(Sys.info()[["sysname"]], Windows = {
      mclapply.hack
    }, Linux = {
      mclapply
    }, Darwin = {
      mclapply
    })

    list.of.seeds.all <- 1:nTrials + Seed

    if (!missing(JobSize))
      JobSize <- JobSize else {
      JobSize <- as.numeric(benchmarkme::get_ram()) * 0.8 * nTrials * (dim(subset(dfx,
        select = -c(Cls)))[1] * dim(subset(dfx, select = -c(Cls)))[2])
    }

    if (nProc > 1) {
      list.of.seeds <- split(list.of.seeds.all, ceiling(seq_along(list.of.seeds.all)/max(nTrials/nProc,
        JobSize)))
    } else {
      list.of.seeds <- split(list.of.seeds.all, 1)
    }

    nlist.of.seeds <- unlist(lapply(list.of.seeds, length))
    ADstatAll <- vector()
    ReducedDataI <- list()
    RemovedDataI <- list()
    for (i in 1:length(list.of.seeds)) {
      ADstat <- vector()
      if (nlist.of.seeds[[i]] * length(Cls) > 6000) {
        ReducedDataMat <- mclapply(1:nlist.of.seeds[i], function(x) {
          set.seed(list.of.seeds[[i]][x])
          sample <- caTools::sample.split(dfx$Cls, SplitRatio = Size/nrow(dfx))
          ReducedDataList <- subset(dfx, sample == TRUE)
          RemovedDataList <- subset(dfx, sample == FALSE)
          ADv <- mapply(CompDistrib, dfx[1:(ncol(dfx) - 1)], ReducedDataList[1:(ncol(ReducedDataList) -
          1)])
          return(list(ReducedDataList = ReducedDataList, RemovedDataList = RemovedDataList,
          ADv = ADv))
        }, mc.cores = nProc)
      } else {
        ReducedDataMat <- lapply(1:nlist.of.seeds[i], function(x) {
          set.seed(list.of.seeds[[i]][x])
          sample <- caTools::sample.split(dfx$Cls, SplitRatio = Size/nrow(dfx))
          ReducedDataList <- subset(dfx, sample == TRUE)
          RemovedDataList <- subset(dfx, sample == FALSE)
          ADv <- mapply(CompDistrib, dfx[1:(ncol(dfx) - 1)], ReducedDataList[1:(ncol(ReducedDataList) -
          1)])
          return(list(ReducedDataList = ReducedDataList, RemovedDataList = RemovedDataList,
          ADv = ADv))
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
  return(list(ReducedData = ReducedData, RemovedData = RemovedData, ReducedInstances = row.names(ReducedData)))
}
