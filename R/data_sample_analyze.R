#' Perform Sampling and Analyze Picked Data Subsets
#'
#' This function performs sampling and analyzes the picked data subsets using the
#' `make_and_analyse_subsample` function.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param TestStat A character string specifying the statistical test to be used for
#'   comparing the distributions.
#' @param Size A numeric value between 0 and 1 representing the desired size of the
#'   downsampled dataset as a proportion of the original dataset.
#' @param list.of.seeds A vector of integer values to be used as the random seeds for
#'   reproducibility.
#' @param PCAimportance A logical value indicating whether to use PCA to identify
#'   relevant variables.
#' @param NonNoiseSelection A logical value indicating whether to use non-uniform
#'   distribution tests to identify relevant variables.
#' @param UniformTestStat A character string specifying the statistical test to be used for
#'   non-uniform variable selection. Available options are: "ks" (Kolmogorov-Smirnov),
#'   "ad" (Anderson-Darling), "kuiper", "cvm" (Cramér-von Mises), "wass" (Wasserstein),
#'   "dts" (Distributional Transform Statistic), "kld" (Kullback-Leibler divergence),
#'   "amrdd" (Average Mean Root of Distributional Differences), and "euc" (Euclidean distance).
#'   Only used when NonNoiseSelection = TRUE (default: "ks").
#' @param UniformThreshold Threshold value for non-uniform variable selection (default: 0.1).
#' @param nProc The number of cores to use for parallel processing.
#' @param CheckRemoved A logical value indicating whether to also optimize the removed part
#'   of the data for distribution equality with the original.
#' @param CheckThreefold A logical value indicating whether to also optimize the reduced part
#'   of the data for distribution equality with the removed part. Ignored when CheckRemoved is FALSE.
#' @param OptimizeBetween A logical value indicating whether to optimize the reduced part
#'   of the data for distribution equality with the removed part. If set, all other comparisions are not performed.
#' @param JobSize Number of seeds to process in each chunk (default: min(50, length(list.of.seeds)/nProc))
#'
#' @return A list of the results from the `make_and_analyse_subsample` function for
#'   each seed in `list.of.seeds`.
#'
#' @importFrom stats prcomp
#' @importFrom pbmcapply pbmclapply
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
sample_and_analyze <- function(DataAndClasses, TestStat, Size, list.of.seeds, PCAimportance, nProc,
                               CheckRemoved, CheckThreefold, OptimizeBetween, JobSize = NULL,
                               NonNoiseSelection = FALSE, UniformTestStat = "ks",
                               UniformThreshold = 0.05) {

  # Set default chunk size if not provided
  if (is.null(JobSize)) {
    JobSize <- min(50, max(1, ceiling(length(list.of.seeds) / max(1, nProc))))
  }

  # Start with all variables (excluding class column)
  all_vars <- names(DataAndClasses)[1:(ncol(DataAndClasses) - 1)]
  selectedVars <- all_vars

  # Variable selection methods
  pca_vars <- character(0)
  nonuniform_vars <- character(0)

  # PCA-based variable selection
  if (PCAimportance && length(list.of.seeds) > 1 && ncol(DataAndClasses) > 2) {
    cat("Applying PCA-based variable selection...\n")
    # Use less memory-intensive PCA options
    pca1 <- prcomp(DataAndClasses[1:(ncol(DataAndClasses) - 1)],
                   retx = FALSE, # Don't store transformed data
                   center = TRUE,
                   scale = TRUE,
                   rank. = min(10, ncol(DataAndClasses) - 1)) # Limit components
    pca_vars <- names(DataAndClasses)[which(names(DataAndClasses) %in% relevant_PCAvariables(pca1))]

    if (length(pca_vars) > 0) {
      selectedVars <- pca_vars
      cat(paste("PCA selected", length(pca_vars), "out of", length(all_vars), "variables\n"))
    } else {
      cat("PCA selection found no variables, using all variables\n")
    }

    rm(pca1) # Clean up PCA object
    gc() # Force garbage collection
  }

  # Non-uniform distribution variable selection
  if (NonNoiseSelection && length(list.of.seeds) > 1 && ncol(DataAndClasses) > 2) {
    cat("Applying non-uniform distribution variable selection...\n")

    nonuniform_vars <- identify_non_noise_variables(
      data_df = DataAndClasses,
      test_stat = UniformTestStat,
      significance_threshold = UniformThreshold,
      verbose = FALSE,
      seed = list.of.seeds[1] # Use first seed for consistent variable selection
    )

    if (length(nonuniform_vars) > 0) {
      # If both PCA and non-uniform selection are used, take the intersection
      if (PCAimportance && length(pca_vars) > 0) {
        combined_vars <- intersect(pca_vars, nonuniform_vars)
        if (length(combined_vars) > 0) {
          selectedVars <- combined_vars
          cat(paste("Combined selection (PCA and Non-uniform) selected", length(combined_vars),
                    "variables\n"))
        } else {
          # If intersection is empty, use union as fallback
          selectedVars <- unique(c(pca_vars, nonuniform_vars))
          cat(paste("No overlap between PCA and non-uniform selection.",
                    "Using union:", length(selectedVars), "variables\n"))
        }
      } else {
        # Only non-uniform selection
        selectedVars <- nonuniform_vars
        cat(paste("Non-uniform selection selected", length(nonuniform_vars),
                  "out of", length(all_vars), "variables\n"))
      }
    } else {
      cat("Non-uniform selection found no variables")
      if (!PCAimportance || length(pca_vars) == 0) {
        cat(", using all variables\n")
      } else {
        cat(", keeping PCA selection\n")
      }
    }
  }

  # Ensure at least one variable is selected
  if (length(selectedVars) == 0) {
    selectedVars <- all_vars[1] # Use first variable as fallback
    message(paste("Warning: No variables selected by any method.",
                  "Using first variable '", selectedVars[1], "' to proceed."))
  }

  # Print final selection summary
  if (PCAimportance || NonNoiseSelection) {
    cat(paste("Final variable selection:", length(selectedVars), "out of",
              length(all_vars), "variables\n"))
    cat(paste("Selected variables:", paste(selectedVars, collapse = ", "), "\n"))
  }

  # Pre-compute data subset for selected variables only
  DataSubset <- DataAndClasses[, c(selectedVars, names(DataAndClasses)[ncol(DataAndClasses)]), drop = FALSE]

  # Memory-optimized sampling function
  make_and_analyse_subsample <- function(DataSubset, TestStat, Size, Seed, selectedVars, CheckRemoved, CheckThreefold, OptimizeBetween) {
    df_reduced <- MakeReducedDataMat(DataSubset, Size, Seed)
    ADv_reduced <- CompareReducedDataMat(DataAndClasses = DataSubset,
                                         ReducedDataList = df_reduced$ReducedDataList,
                                         TestStat = TestStat)

    if (CheckRemoved || OptimizeBetween) {
      ADv_removed <- CompareReducedDataMat(DataAndClasses = DataSubset,
                                           ReducedDataList = df_reduced$RemovedDataList,
                                           TestStat = TestStat)
    } else {
      ADv_removed <- rep(NA_real_, length(selectedVars))
      names(ADv_removed) <- selectedVars
    }

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

  # Process in chunks to reduce memory pressure
  process_chunk <- function(seed_chunk, DataSubset, TestStat, Size, selectedVars, CheckRemoved, CheckThreefold, OptimizeBetween, use_parallel = FALSE, cores = 1) {
    if (use_parallel && cores > 1 && length(seed_chunk) > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        # Use foreach for Windows
        doParallel::registerDoParallel(min(cores, length(seed_chunk)))
        seed <- integer()
        result <- foreach::foreach(seed = seed_chunk,
                                   .combine = 'c',
                                   .maxcombine = length(seed_chunk),
                                   .multicombine = TRUE,
                                   .packages = c()) %dopar% {
          list(make_and_analyse_subsample(DataSubset, TestStat, Size, seed, selectedVars, CheckRemoved, CheckThreefold, OptimizeBetween))
        }
        doParallel::stopImplicitCluster()
      } else {
        # Use mclapply for Unix-like systems
        result <- pbmcapply::pbmclapply(
          seed_chunk,
          function(seed) make_and_analyse_subsample(DataSubset, TestStat, Size, seed, selectedVars, CheckRemoved, CheckThreefold, OptimizeBetween),
          mc.cores = min(cores, length(seed_chunk)),
          mc.preschedule = TRUE # Better load balancing
        )
      }
    } else {
      # Sequential processing
      result <- lapply_with_bar(
        seed_chunk,
        function(seed) make_and_analyse_subsample(DataSubset, TestStat, Size, seed, selectedVars, CheckRemoved, CheckThreefold, OptimizeBetween)
      )
    }
    return(result)
  }

  # Split seeds into chunks
  seed_chunks <- split(list.of.seeds, ceiling(seq_along(list.of.seeds) / JobSize))

  # Initialize result list
  ReducedDataMat <- vector("list", length(list.of.seeds))
  current_index <- 1

  # Process chunks
  for (chunk in seed_chunks) {
    chunk_result <- process_chunk(
      seed_chunk = chunk,
      DataSubset = DataSubset,
      TestStat = TestStat,
      Size = Size,
      selectedVars = selectedVars,
      CheckRemoved = CheckRemoved,
      CheckThreefold = CheckThreefold,
      OptimizeBetween = OptimizeBetween,
      use_parallel = nProc > 1,
      cores = nProc
    )

    # Store results
    end_index <- current_index + length(chunk) - 1
    ReducedDataMat[current_index:end_index] <- chunk_result
    current_index <- end_index + 1

    # Force garbage collection after each chunk
    gc()
  }

  return(ReducedDataMat)
}