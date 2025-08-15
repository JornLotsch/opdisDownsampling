#' Process Seeds in Chunks with Memory Management
#'
#' This function processes seeds in chunks to manage memory usage effectively
#' during parallel processing of downsampling trials.
#'
#' @param list.of.seeds A vector of integer values to be used as the random seeds
#'   for reproducibility.
#' @param DataSubset A data frame containing the pre-filtered data subset with
#'   selected variables and class column.
#' @param processing_params A list containing parameters for the analysis function:
#'   \itemize{
#'     \item TestStat: Statistical test for comparing distributions
#'     \item Size: Desired size of downsampled dataset
#'     \item selectedVars: Vector of selected variable names
#'     \item CheckRemoved: Logical for optimizing removed data
#'     \item CheckThreefold: Logical for three-way optimization
#'     \item OptimizeBetween: Logical for between-group optimization
#'   }
#' @param JobSize Number of seeds to process in each chunk.
#' @param nProc Number of processor cores to use for parallel processing.
#'
#' @return A list of results from all seeds processed across all chunks.
#'
#' @details This function manages memory by:
#' \itemize{
#'   \item Splitting seeds into manageable chunks
#'   \item Processing each chunk separately
#'   \item Forcing garbage collection after each chunk
#'   \item Maintaining result ordering across chunks
#' }
#'
#' @keywords internal
process_seeds_in_chunks <- function(list.of.seeds, DataSubset, processing_params, JobSize, nProc) {

  # Split seeds into chunks
  seed_chunks <- split(list.of.seeds, ceiling(seq_along(list.of.seeds) / JobSize))

  # Initialize result list
  ReducedDataMat <- vector("list", length(list.of.seeds))
  current_index <- 1

  # Process chunks
  for (chunk in seed_chunks) {
    chunk_result <- process_single_chunk(chunk, DataSubset, processing_params, nProc)

    # Store results
    end_index <- current_index + length(chunk) - 1
    ReducedDataMat[current_index:end_index] <- chunk_result
    current_index <- end_index + 1

    # Force garbage collection after each chunk
    gc()
  }

  return(ReducedDataMat)
}

#' Process a Single Chunk of Seeds
#'
#' This function processes a single chunk of seeds, choosing between parallel
#' and sequential processing based on the number of available processors and
#' chunk size.
#'
#' @param seed_chunk A vector of seeds for this specific chunk.
#' @param DataSubset A data frame containing the pre-filtered data subset.
#' @param processing_params A list containing parameters for the analysis function.
#' @param nProc Number of processor cores available for parallel processing.
#'
#' @return A list of results for this chunk of seeds.
#'
#' @details The function automatically selects the processing method:
#' \itemize{
#'   \item Parallel processing when nProc > 1 and chunk has multiple seeds
#'   \item Sequential processing otherwise
#' }
#'
#' @keywords internal
process_single_chunk <- function(seed_chunk, DataSubset, processing_params, nProc) {

  use_parallel <- nProc > 1 && length(seed_chunk) > 1

  if (use_parallel) {
    return(process_chunk_parallel(seed_chunk, DataSubset, processing_params, nProc))
  } else {
    return(process_chunk_sequential(seed_chunk, DataSubset, processing_params))
  }
}

#' Process Chunk Using Parallel Processing
#'
#' This function processes a chunk of seeds using parallel processing, with
#' platform-specific optimization for Windows and Unix-like systems.
#'
#' @param seed_chunk A vector of seeds to process in parallel.
#' @param DataSubset A data frame containing the pre-filtered data subset.
#' @param processing_params A list containing parameters for make_and_analyse_subsample.
#' @param nProc Number of processor cores to use.
#'
#' @return A list of results from parallel processing.
#'
#' @details Platform-specific implementations:
#' \itemize{
#'   \item Windows: Uses foreach with doParallel backend
#'   \item Unix-like: Uses pbmclapply with multicore backend
#' }
#'
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom pbmcapply pbmclapply
#'
#' @keywords internal
process_chunk_parallel <- function(seed_chunk, DataSubset, processing_params, nProc) {

  if (Sys.info()[["sysname"]] == "Windows") {
    # Use foreach for Windows
    doParallel::registerDoParallel(min(nProc, length(seed_chunk)))

    # Initialize variable for foreach
    seed <- integer()

    result <- foreach::foreach(
      seed = seed_chunk,
      .combine = 'c',
      .maxcombine = length(seed_chunk),
      .multicombine = TRUE,
      .packages = c()
    ) %dopar% {
      list(make_and_analyse_subsample(
        DataSubset = DataSubset,
        TestStat = processing_params$TestStat,
        Size = processing_params$Size,
        Seed = seed,
        selectedVars = processing_params$selectedVars,
        CheckRemoved = processing_params$CheckRemoved,
        CheckThreefold = processing_params$CheckThreefold,
        OptimizeBetween = processing_params$OptimizeBetween
      ))
    }

    doParallel::stopImplicitCluster()

  } else {
    # Use mclapply for Unix-like systems
    result <- pbmcapply::pbmclapply(
      seed_chunk,
      function(seed) {
        make_and_analyse_subsample(
          DataSubset = DataSubset,
          TestStat = processing_params$TestStat,
          Size = processing_params$Size,
          Seed = seed,
          selectedVars = processing_params$selectedVars,
          CheckRemoved = processing_params$CheckRemoved,
          CheckThreefold = processing_params$CheckThreefold,
          OptimizeBetween = processing_params$OptimizeBetween
        )
      },
      mc.cores = min(nProc, length(seed_chunk)),
      mc.preschedule = TRUE  # Better load balancing
    )
  }

  return(result)
}

#' Process Chunk Using Sequential Processing
#'
#' This function processes a chunk of seeds sequentially with a progress bar
#' for user feedback during long-running operations.
#'
#' @param seed_chunk A vector of seeds to process sequentially.
#' @param DataSubset A data frame containing the pre-filtered data subset.
#' @param processing_params A list containing parameters for make_and_analyse_subsample.
#'
#' @return A list of results from sequential processing.
#'
#' @details Uses lapply_with_bar for progress tracking during sequential
#' processing of seeds.
#'
#' @keywords internal
process_chunk_sequential <- function(seed_chunk, DataSubset, processing_params) {

  result <- lapply_with_bar(
    seed_chunk,
    function(seed) {
      make_and_analyse_subsample(
        DataSubset = DataSubset,
        TestStat = processing_params$TestStat,
        Size = processing_params$Size,
        Seed = seed,
        selectedVars = processing_params$selectedVars,
        CheckRemoved = processing_params$CheckRemoved,
        CheckThreefold = processing_params$CheckThreefold,
        OptimizeBetween = processing_params$OptimizeBetween
      )
    }
  )

  return(result)
}