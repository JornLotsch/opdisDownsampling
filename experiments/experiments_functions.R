# ===============================================================================
# FUNCTIONS USED IN THE ARTIFICIAL DATA SETS GENERATION AND EVALUATION SCRIPT
# ===============================================================================

# Attempt to set working directory to script location (for RStudio)
tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable("getSourceEditorContext")) {
    script_path <- rstudioapi::getSourceEditorContext()$path
    if (!is.null(script_path) && nzchar(script_path)) {
      setwd(dirname(script_path))
      cat("Working directory set to script location:", getwd(), "\n")
    }
  }
}, error = function(e) {
  message("Unable to set working directory automatically: ", e$message)
})

# Change directory if it exists and is different from current
if (dir.exists(EXPERIMENTS_DIR)) {
  tryCatch({
    setwd(EXPERIMENTS_DIR)
    cat("Working directory changed to:", getwd(), "\n")
  }, error = function(e) {
    warning("Failed to change to experiments directory: ", e$message)
    cat("Continuing with current directory:", getwd(), "\n")
  })
} else if (!dir.exists(EXPERIMENTS_DIR)) {
  warning("Experiments directory does not exist: ", EXPERIMENTS_DIR)
  cat("Continuing with current directory:", getwd(), "\n")
}

# Verify current working directory
cat("Current working directory:", getwd(), "\n")


# ===============================================================================
# LOAD REQUIRED LIBRARIES
# ===============================================================================
library(stats) # Statistical functions
library(ggplot2) # Data visualization
library(ggthemes) # Additional ggplot2 themes
library(dplyr) # Data manipulation
library(parallel) # Parallel processing
library(pbmcapply) # Progress bars for parallel processing
library(cowplot) # Plot arrangement
library(reshape2) # Data reshaping
library(tidyr) # Data tidying
library(ABCanalysis) # Data categorization

# ============================================================================
# UTILITY FUNCTIONS
# ===============================================================================

#' Perform Data Splitting for Single Iteration
#'
#' This function performs a single iteration of data splitting using the specified method.
#' Designed to be used within lapply loops.
#' @importFrom base %in%
#' @param data_list List containing 'Data' and 'Cls' components
#' @param downsampling_size Target size for split data
#' @param iteration Current iteration number (for seed generation)
#' @param split_method Character string specifying splitting method ("opdis", "random", etc.)
#' @param split_params List of method-specific parameters
#'
#' @return Result from the splitting method (structure depends on method)
#'
perform_single_split <- function(data_list, downsampling_size, iteration,
                                 split_method = "opdis", split_params = list()) {

  switch(split_method,
         "opdis" = {
    # Set default parameters for OPDIS - complete parameter list
    default_params <- list(
             nTrials = 1,
             CheckRemoved = FALSE,
             CheckThreefold = FALSE,
             OptimizeBetween = FALSE,
             NonNoiseSelection = FALSE,
             TestStat = "ad",
             WorstSample = FALSE,
             PCAimportance = FALSE,
             UniformTestStat = "ks",
             UniformThreshold = 0.05,
             JobSize = 0,
             verbose = FALSE
           )

    # Merge with provided parameters
    params <- modifyList(default_params, split_params)

    # Generate seed exactly like original code
    seed <- iteration + (iteration - 1) * 1e6

    # Calculate MaxCores exactly like original code
    max_cores <- min(params$nTrials, parallel::detectCores() - 1)

    # Call opdisDownsampling with complete parameters
    opdisDownsampling::opdisDownsampling(
             Data = data_list$Data,
             Cls = data_list$Cls,
             Size = downsampling_size,
             Seed = seed,
             nTrials = params$nTrials,
             CheckRemoved = params$CheckRemoved,
             CheckThreefold = params$CheckThreefold,
             OptimizeBetween = params$OptimizeBetween,
             MaxCores = max_cores,
             TestStat = params$TestStat,
             PCAimportance = params$PCAimportance,
             NonNoiseSelection = params$NonNoiseSelection,
             WorstSample = params$WorstSample,
             UniformTestStat = params$UniformTestStat,
             UniformThreshold = params$UniformThreshold,
             JobSize = params$JobSize,
             verbose = params$verbose
           )
  },

         "random" = {
    # Random splitting implementation
    default_params <- list(stratified = TRUE)
    params <- modifyList(default_params, split_params)

    # Set seed for reproducibility
    set.seed(iteration + (iteration - 1) * 1e6)

    n <- nrow(data_list$Data)
    target_size <- if (downsampling_size < 1) {
      round(downsampling_size * n)
    } else {
      as.integer(downsampling_size)
    }

    if (params$stratified && length(unique(data_list$Cls)) > 1) {
      # Stratified sampling
      reduced_indices <- unlist(lapply(unique(data_list$Cls), function(cls) {
        cls_indices <- which(data_list$Cls == cls)
        cls_target <- round(length(cls_indices) * target_size / n)
        if (cls_target > 0) {
          sample(cls_indices, min(cls_target, length(cls_indices)))
        } else {
          integer(0)
        }
      }))
    } else {
      # Simple random sampling
      reduced_indices <- sample(n, target_size)
    }

    removed_indices <- setdiff(1:n, reduced_indices)

    # Return in opdis format for consistency
    list(
             ReducedData = data.frame(
               data_list$Data[reduced_indices,, drop = FALSE],
               Cls = data_list$Cls[reduced_indices]
             ),
             RemovedData = data.frame(
               data_list$Data[removed_indices,, drop = FALSE],
               Cls = data_list$Cls[removed_indices]
             )
           )
  },

         stop("Unsupported split method: ", split_method)
  )
}

#' Validate input data structure and return feature columns
#' @param data_df Data frame to validate
#' @param class_name Name of the class column
#' @return Character vector of feature column names
validate_input_data <- function(data_df, class_name) {
  if (!is.data.frame(data_df)) {
    stop("data_df must be a data frame")
  }

  if (!class_name %in% names(data_df)) {
    stop(paste("Class column", class_name, "not found in data"))
  }

  if (nrow(data_df) == 0) {
    stop("data_df is empty")
  }

  feature_cols <- setdiff(names(data_df), class_name)
  if (length(feature_cols) == 0) {
    stop("No feature columns found")
  }

  # Check if feature columns are numeric
  numeric_features <- sapply(data_df[feature_cols], is.numeric)
  if (!all(numeric_features)) {
    non_numeric <- feature_cols[!numeric_features]
    stop(paste("Non-numeric feature columns found:", paste(non_numeric, collapse = ", ")))
  }

  return(feature_cols)
}

#' Check for redundancy in sampling results
#' @param downsampled_results List of downsampling results
#' @param overlap_threshold Threshold for considering high overlap (default: 0.95)
#' @param check_identical Whether to check for identical datasets (default: TRUE)
#' @param sample_size For large datasets, sample this many rows for comparison (default: NULL for all)
#' @return Logical indicating if redundancy was found
check_redundancy <- function(downsampled_results, overlap_threshold = 0.95,
                             check_identical = TRUE, sample_size = NULL) {
  cat("Checking for redundancy in sampling results...\n")

  # Extract reduced data from all iterations
  reduced_data_list <- lapply(downsampled_results, function(res) {
    if (is.null(res$reduced_data)) return(NULL)
    # Convert to matrix for comparison, excluding class column
    data_cols <- setdiff(names(res$reduced_data), "Cls")
    as.matrix(res$reduced_data[, data_cols, drop = FALSE])
  })

  # Remove NULL entries
  reduced_data_list <- reduced_data_list[!sapply(reduced_data_list, is.null)]

  if (length(reduced_data_list) <= 1) {
    cat("Insufficient data for redundancy check.\n")
    return(FALSE)
  }

  n_iterations <- length(reduced_data_list)
  redundancy_found <- FALSE

  # Early termination if datasets have different dimensions
  dimensions <- sapply(reduced_data_list, dim)
  unique_dims <- unique(t(dimensions))

  if (nrow(unique_dims) > 1) {
    cat("Datasets have different dimensions - no redundancy possible.\n")
    return(FALSE)
  }

  # Optimization 1: Use hashing for identical dataset detection
  if (check_identical) {
    cat("Checking for identical datasets using hash comparison...\n")

    # Calculate hash for each dataset for quick identical comparison
    dataset_hashes <- sapply(reduced_data_list, function(data) {
      if (is.null(sample_size) || nrow(data) <= sample_size) {
        digest::digest(data, algo = "md5")
      } else {
        # For very large datasets, hash a sample
        sample_rows <- sample(nrow(data), sample_size)
        digest::digest(data[sample_rows,, drop = FALSE], algo = "md5")
      }
    })

    # Find duplicate hashes
    hash_duplicates <- duplicated(dataset_hashes)
    if (any(hash_duplicates)) {
      duplicate_pairs <- which(hash_duplicates)
      for (dup_idx in duplicate_pairs) {
        original_idx <- which(dataset_hashes == dataset_hashes[dup_idx])[1]
        cat(paste("Warning: Identical samples found between iterations", original_idx, "and", dup_idx, "\n"))
        redundancy_found <- TRUE
      }
    }
  }

  # Optimization 2: Vectorized overlap calculation for datasets of same size
  if (length(unique_dims) == 1) {
    cat("Checking for high overlap using vectorized comparison...\n")

    # Convert list to 3D array for vectorized operations
    n_rows <- unique_dims[1, 1]
    n_cols <- unique_dims[1, 2]

    # Only proceed if datasets are not too large for memory
    total_memory_needed <- n_iterations * n_rows * n_cols * 8 # 8 bytes per double
    memory_limit <- 1e9 # 1GB limit

    if (total_memory_needed < memory_limit) {
      # Stack all datasets into a 3D array
      data_array <- array(dim = c(n_rows, n_cols, n_iterations))
      for (i in seq_along(reduced_data_list)) {
        data_array[,, i] <- reduced_data_list[[i]]
      }

      # Vectorized pairwise comparison
      for (i in 1:(n_iterations - 1)) {
        # Calculate overlaps with all subsequent datasets at once
        remaining_indices <- (i + 1):n_iterations

        # Vectorized overlap calculation
        overlaps <- sapply(remaining_indices, function(j) {
          # Count matching rows between datasets i and j
          matches <- sum(apply(data_array[,, i], 1, function(row1) {
            any(apply(data_array[,, j], 1, function(row2) {
              all(abs(row1 - row2) < 1e-10)
            }))
          }))
          matches / n_rows
        })

        # Check for high overlaps
        high_overlap_indices <- which(overlaps > overlap_threshold)
        if (length(high_overlap_indices) > 0) {
          for (idx in high_overlap_indices) {
            j <- remaining_indices[idx]
            cat(paste("Warning: High overlap (", round(overlaps[idx] * 100, 1),
                      "%) found between iterations", i, "and", j, "\n"))
            redundancy_found <- TRUE
          }
        }
      }
    } else {
      cat("Datasets too large for vectorized comparison, using optimized pairwise approach...\n")

      # Optimization 3: Sample-based comparison for very large datasets
      comparison_sample_size <- min(sample_size %||% 1000, n_rows)

      for (i in 1:(n_iterations - 1)) {
        for (j in (i + 1):n_iterations) {
          # Sample rows for comparison if datasets are large
          if (n_rows > comparison_sample_size) {
            sample_rows <- sample(n_rows, comparison_sample_size)
            data_i_sample <- reduced_data_list[[i]][sample_rows,, drop = FALSE]
            data_j_sample <- reduced_data_list[[j]][sample_rows,, drop = FALSE]
          } else {
            data_i_sample <- reduced_data_list[[i]]
            data_j_sample <- reduced_data_list[[j]]
          }

          # Fast overlap calculation using sample
          overlap <- sum(apply(data_i_sample, 1, function(row1) {
            any(apply(data_j_sample, 1, function(row2) {
              all(abs(row1 - row2) < 1e-10)
            }))
          })) / nrow(data_i_sample)

          if (overlap > overlap_threshold) {
            cat(paste("Warning: High overlap (", round(overlap * 100, 1),
                      "%) found between iterations", i, "and", j,
                      ifelse(n_rows > comparison_sample_size, " (based on sample)", ""), "\n"))
            redundancy_found <- TRUE
          }
        }
      }
    }
  }

  if (!redundancy_found) {
    cat("No significant redundancy found across sampling iterations.\n")
  }

  return(redundancy_found)
}

# Helper function for null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Error handling wrapper for statistical tests
#' @param test_func Function to execute
#' @param error_context Description of the context for error messages
#' @param default_value Default value to return on error
safe_test_execution <- function(test_func, error_context, default_value = NA) {
  tryCatch({
    test_func()
  }, error = function(e) {
    warning(paste(error_context, "failed:", e$message), call. = FALSE)
    default_value
  })
}

# ===============================================================================
# MAIN EXPERIMENT FUNCTION
# ===============================================================================

#' Run downsampling experiment with provided data
#' @param data_df Data frame containing the dataset to analyze
#' @param class_name Name of the class column
#' @param CheckRemoved Logical; optimize for removed parts
#' @param CheckThreefold Logical; optimize for removed and between reduced and removed parts
#' @param OptimizeBetween Logical; optimize between reduced and removed parts
#' @param NonNoiseSelection Logical; use non-uniform variable selection
#' @return List containing comprehensive analysis results
run_experiment <- function(data_df, class_name, CheckRemoved = FALSE, CheckThreefold = FALSE,
                           OptimizeBetween = FALSE, NonNoiseSelection = FALSE,
                           PCAimportance = FALSE, UniformTestStat = "ks",
                           UniformThreshold = 0.05, JobSize = 0, verbose = FALSE,
                           WorstSample = FALSE) {
  # Input validation
  feature_cols <- validate_input_data(data_df, class_name)

  # Convert data frame to list format expected by opdisDownsampling
  data_list <- list(
    Data = data_df[, feature_cols, drop = FALSE],
    Cls = data_df[[class_name]]
  )

  cat("Starting experiment with provided data...\n")
  cat(paste("Data dimensions:", nrow(data_list$Data), "samples x", ncol(data_list$Data), "features\n"))
  cat(paste("Classes:", paste(unique(data_list$Cls), collapse = ", "), "\n"))
  cat(paste("Class distribution:", paste(table(data_list$Cls), collapse = ", "), "\n"))

  # Print current parameter combination for tracking progress
  cat(paste0("Parameters: nTrials = ", nTrials,
            ", CheckRemoved = ", CheckRemoved,
            ", CheckThreefold = ", CheckThreefold,
            ", OptimizeBetween = ", OptimizeBetween,
            ", NonNoiseSelection = ", NonNoiseSelection, "\n"))

  #' Compute U-test p-values for each variable distinguishing classes
  #' Uses Wilcoxon rank-sum test to assess class discrimination ability
  calculate_Utest_pvalues <- function(Data, Cls) {
    # Remove class column if it exists in Data
    if ("Cls" %in% names(Data)) {
      Data <- Data[, !names(Data) %in% "Cls", drop = FALSE]
    }

    # Ensure Data is a data frame or matrix
    if (is.vector(Data)) {
      Data <- as.data.frame(Data)
    }

    # Calculate p-values for each variable
    p_values <- sapply(1:ncol(Data), function(i) {
      var_data <- Data[, i]
      class1_data <- var_data[Cls == 1]
      class2_data <- var_data[Cls == 2]

      # Remove NA values
      class1_data <- class1_data[!is.na(class1_data)]
      class2_data <- class2_data[!is.na(class2_data)]

      # Check if we have enough data
      if (length(class1_data) < 2 || length(class2_data) < 2) {
        return(1.0) # Return non-significant p-value if insufficient data
      }

      # Perform Mann-Whitney U test (Wilcoxon rank-sum test)
      tryCatch({
        test_result <- wilcox.test(class1_data, class2_data, exact = FALSE)
        return(test_result$p.value)
      }, error = function(e) {
        return(1.0) # Return non-significant p-value on error
      })
    })

    # Return structured result
    return(list(
      Variable_No = 1:ncol(Data),
      p_values = p_values
    ))
  }

  calculate_regression_pvalues <- function(Data, Cls) {
    # Remove class column if it exists in Data
    if ("Cls" %in% names(Data)) {
      Data <- Data[, !names(Data) %in% "Cls", drop = FALSE]
    }

    # Ensure Data is a data frame or matrix
    if (is.vector(Data)) {
      Data <- as.data.frame(Data)
    }

    # Calculate p-values for each variable using logistic regression
    p_values <- sapply(1:ncol(Data), function(i) {
      var_data <- Data[, i]

      # Remove NA values
      valid_indices <- !is.na(var_data) & !is.na(Cls)
      if (sum(valid_indices) < 4) {
        # Need at least 4 observations for regression
        return(1.0)
      }

      var_clean <- var_data[valid_indices]
      cls_clean <- Cls[valid_indices]

      # Check if we have both classes
      if (length(unique(cls_clean)) < 2) {
        return(1.0)
      }

      # Perform logistic regression
      tryCatch({
        # Convert to binary (0/1) for logistic regression
        binary_cls <- as.numeric(cls_clean == max(cls_clean))
        model <- glm(binary_cls ~ var_clean, family = binomial())
        summary_model <- summary(model)

        # Extract p-value for the variable coefficient (second coefficient)
        if (nrow(summary_model$coefficients) >= 2) {
          return(summary_model$coefficients[2, 4]) # p-value is in column 4
        } else {
          return(1.0)
        }
      }, error = function(e) {
        return(1.0) # Return non-significant p-value on error
      })
    })

    # Return structured result
    return(list(
      Variable_No = 1:ncol(Data),
      p_values = p_values
    ))
  }


  #' Plot significance comparison across multiple downsampling iterations
  #' Creates faceted visualization comparing original, reduced, and removed datasets

  plot_significance_comparison_multi <- function(original_pvals, reduced_pvals_list, removed_pvals_list) {
    # Data preparation with better error handling
    prepare_plot_data <- function(pvals, dataset_name, iteration_prefix = "Iter_") {
      if (is.data.frame(pvals)) {
        # Handle data frame format (like pval_results)
        # Check if it has rownames as Variable_No or a Variable_No column
        if ("Variable_No" %in% colnames(pvals)) {
          var_no <- pvals$Variable_No
          p_values <- pvals$p_values
        } else {
          # Use rownames as Variable_No
          var_no <- as.numeric(gsub("\\D", "", rownames(pvals))) # Extract numbers from rownames
          if (all(is.na(var_no))) {
            var_no <- seq_len(nrow(pvals)) # Fallback to sequential numbering
          }
          p_values <- pvals$p_values
        }

        data.frame(
          Variable_No = var_no,
          neg_log_p = -log10(pmax(p_values, 1e-100)), # Avoid -Inf
          Dataset = dataset_name,
          Iteration = dataset_name
        )
      } else if (is.list(pvals) && !is.null(pvals$p_values)) {
        # Single dataset with list structure (has $Variable_No and $p_values)
        data.frame(
          Variable_No = pvals$Variable_No,
          neg_log_p = -log10(pmax(pvals$p_values, 1e-100)), # Avoid -Inf
          Dataset = dataset_name,
          Iteration = dataset_name
        )
      } else if (is.list(pvals) && is.null(pvals$p_values)) {
        # Multiple iterations - this is a list of results
        result_list <- lapply(seq_along(pvals), function(i) {
          if (is.null(pvals[[i]]) || is.null(pvals[[i]]$p_values)) return(NULL)
          data.frame(
            Variable_No = pvals[[i]]$Variable_No,
            neg_log_p = -log10(pmax(pvals[[i]]$p_values, 1e-100)), # Avoid -Inf
            Dataset = dataset_name,
            Iteration = paste0(iteration_prefix, i)
          )
        })

        # Remove NULL results and combine
        result_list <- result_list[!sapply(result_list, is.null)]
        if (length(result_list) > 0) {
          do.call(rbind, result_list)
        } else {
          NULL
        }
      } else {
        # Handle unexpected format
        warning(paste("Unexpected format for", dataset_name, "p-values"))
        return(NULL)
      }
    }

    # Prepare data for each dataset type
    original_data <- prepare_plot_data(pvals = original_pvals, dataset_name = "Original", iteration_prefix = "Original")
    reduced_data <- prepare_plot_data(reduced_pvals_list, "Reduced")
    removed_data <- prepare_plot_data(removed_pvals_list, "Removed")

    # Remove NULL datasets
    all_datasets <- list(original_data, reduced_data, removed_data)
    all_datasets <- all_datasets[!sapply(all_datasets, is.null)]

    if (length(all_datasets) == 0) {
      stop("No valid datasets to plot")
    }

    # Combine all data and set factor levels for consistent plotting
    combined_data <- do.call(rbind, all_datasets)
    combined_data$Dataset <- factor(combined_data$Dataset, levels = c("Original", "Reduced", "Removed"))

    # Define plot theme for consistent styling
    plot_theme <- theme_light() +
      theme(
        legend.position = "bottom",
        legend.background = element_rect(fill = alpha("white", 0.5)),
        strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black"),
        plot.subtitle = element_text(size = 9),
        panel.grid.minor = element_blank()
      )

    # Create the multi-panel
    ggplot(combined_data, aes(x = Variable_No, y = neg_log_p, shape = as.factor(Iteration))) +
      geom_point(data = subset(combined_data, Dataset == "Original"),
                 color = "steelblue", alpha = 0.8, size = 1.5) +
      geom_point(data = subset(combined_data, Dataset == "Reduced"),
                 aes(color = Dataset), alpha = 0.6, size = 1.2) +
      geom_point(data = subset(combined_data, Dataset == "Removed"),
                 aes(color = Dataset), alpha = 0.6, size = 1.2) +
      geom_hline(yintercept = -log10(significance_threshold),
                 color = "salmon", linetype = "dashed", linewidth = 1) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 3) +
      scale_color_colorblind() +
      plot_theme +
      labs(
        x = "Variable Number",
        y = ifelse(use_regression_for_p, "-log10(p-value) regression", "-log10(p-value) U-test"),
        title = paste0("Statistical significance comparison: Original vs ",
                       length(reduced_pvals_list), " Downsampling Iterations"),
        subtitle = paste0("Based on ", format(nTrials, scientific = FALSE),
                          " trials per iteration; CheckRemoved = ", CheckRemoved,
                          ", CheckThreefold = ", CheckThreefold, ", OptimizeBetween = ", OptimizeBetween,
                          ", NonNoiseSelection = ", NonNoiseSelection),
        color = "Dataset Type"
      ) +
      guides(shape = "none", color = "none")
  }

  # ---------------------------------------------------------------------------
  # MAIN EXPERIMENT EXECUTION
  # ---------------------------------------------------------------------------
  cat("Calculating p-values for original dataset...\n")
  # Calculate statistical significance for all variables in original dataset
  pval_results <- if (use_regression_for_p) {
    calculate_regression_pvalues(data_list$Data, data_list$Cls)
  } else {
    calculate_Utest_pvalues(data_list$Data, data_list$Cls)
  }

  # Validate p-value results
  if (any(is.na(pval_results$p_values))) {
    n_na <- sum(is.na(pval_results$p_values))
    cat(paste("Warning:", n_na, "out of", nrow(pval_results), "p-values are NA\n"))
  }

  cat(paste0("Performing ", nSamples, " downsampling iterations...\n"))
  # Perform multiple downsampling iterations with different seeds
  downsampled_results <- lapply(1:nSamples, function(i) {
    cat(paste0("  Iteration ", i, "/", nSamples, "\n"))

    result <- tryCatch({
      perform_single_split(
        data_list = data_list,
        downsampling_size = downsampling_size,
        iteration = i,
        split_method = split_method,
        split_params = list(
          nTrials = nTrials,
          CheckRemoved = CheckRemoved,
          CheckThreefold = CheckThreefold,
          OptimizeBetween = OptimizeBetween,
          NonNoiseSelection = NonNoiseSelection,
          TestStat = TestStat,
          WorstSample = WorstSample,
          PCAimportance = PCAimportance, # Add this if it's available in scope
      # Add the missing parameters:
          UniformTestStat = UniformTestStat, # Default: "ks"
          UniformThreshold = UniformThreshold, # Default: 0.05
          JobSize = JobSize, # Default: 0
          verbose = verbose # Default: FALSE
        )
      )
    }, error = function(e) {
      warning(paste("Error in iteration", i, ":", e$message), call. = FALSE)
      return(NULL)
    })

    if (is.null(result)) return(NULL)

    # Return both reduced (kept) and removed data with class labels
    list(
      reduced_data = cbind(Cls = result$ReducedData$Cls,
                           result$ReducedData[, - ncol(result$ReducedData), drop = FALSE]),
      removed_data = cbind(Cls = result$RemovedData$Cls,
                           result$RemovedData[, - ncol(result$RemovedData), drop = FALSE])
    )
  })

  # Remove NULL results from failed iterations
  downsampled_results <- downsampled_results[!sapply(downsampled_results, is.null)]
  names(downsampled_results) <- paste0("Iteration", seq_along(downsampled_results))

  if (length(downsampled_results) == 0) {
    stop("All downsampling iterations failed")
  }

  if (length(downsampled_results) < nSamples) {
    cat(paste("Warning: Only", length(downsampled_results), "out of", nSamples, "iterations succeeded\n"))
  }

  # ---------------------------------------------------------------------------
  # REDUNDANCY CHECK
  # ---------------------------------------------------------------------------
  redundancy_found <- check_redundancy(downsampled_results)

  # ---------------------------------------------------------------------------
  # Calculating p-values in sampled data
  # ---------------------------------------------------------------------------

  cat("Calculating p-values for reduced datasets...\n")
  # Calculate p-values for each reduced dataset
  pval_reduced_list <- lapply(downsampled_results, function(res) {
    if (is.null(res$reduced_data)) return(NULL)
    if (use_regression_for_p) {
      calculate_regression_pvalues(res$reduced_data, res$reduced_data$Cls)
    } else {
      calculate_Utest_pvalues(res$reduced_data, res$reduced_data$Cls)
    }
  })

  cat("Calculating p-values for removed datasets...\n")
  # Calculate p-values for each removed dataset
  pval_removed_list <- lapply(downsampled_results, function(res) {
    if (is.null(res$removed_data)) return(NULL)
    if (use_regression_for_p) {
      calculate_regression_pvalues(res$removed_data, res$removed_data$Cls)
    } else {
      calculate_Utest_pvalues(res$removed_data, res$removed_data$Cls)
    }
  })

  # Remove NULL results
  pval_reduced_list <- pval_reduced_list[!sapply(pval_reduced_list, is.null)]
  pval_removed_list <- pval_removed_list[!sapply(pval_removed_list, is.null)]

  cat("Creating visualization...\n")
  # Create main significance comparison visualization
  comparison_plot <- plot_significance_comparison_multi(original_pvals = pval_results,
                                                        reduced_pvals_list = pval_reduced_list,
                                                        removed_pvals_list = pval_removed_list)

  # ---------------------------------------------------------------------------
  # JACCARD INDEX ANALYSIS (Feature Selection Overlap)
  # ---------------------------------------------------------------------------

  # Calculate Jaccard index for feature selection overlap between sets
  calculate_jaccard_index <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    if (union == 0) return(1) # Both sets are empty - perfect agreement
    return(intersection / union)
  }

  # P-value calculation with Jaccard analysis
  # Compares feature selection overlap between: Original-Reduced, Original-Removed, Reduced-Removed
  calculate_pvals_with_jaccard_threefold <- function(pval_results, pval_reduced_list, pval_removed_list,
                                                     significance_threshold = 0.05) {

    # Apply FDR correction for multiple testing
    fdr_adjusted_pvals <- p.adjust(pval_results$p_values, method = "fdr")

    cat("Calculating Jaccard indices for feature selection overlap...\n")

    # Identify significant features for original data using both raw and FDR-corrected p-values
    if (!use_ABC_for_feature_selection || !requireNamespace("ABCanalysis", quietly = TRUE)) {
      original_significant_raw <- which(pval_results$p_values < significance_threshold)
      original_significant_fdr <- which(fdr_adjusted_pvals < significance_threshold)
    } else {
      # Use ABCanalysis for feature selection based on -log10 transformed p-values
      original_significant_raw <- safe_test_execution(
        function() ABCanalysis::ABCanalysis(-log10(pmax(pval_results$p_values, 1e-100)))$Aind,
        "ABCanalysis for raw p-values",
        which(pval_results$p_values < significance_threshold)
      )

      original_significant_fdr <- safe_test_execution(
        function() ABCanalysis::ABCanalysis(-log10(pmax(fdr_adjusted_pvals, 1e-100)))$Aind,
        "ABCanalysis for FDR-corrected p-values",
        which(fdr_adjusted_pvals < significance_threshold)
      )
    }

    # Initialize results data frame for all three comparison types
    jaccard_results <- data.frame(
      Iteration = integer(),
      Comparison = character(),
      Jaccard_Raw = numeric(),
      Jaccard_FDR = numeric(),
      stringsAsFactors = FALSE
    )

    # Process each downsampling iteration
    for (i in seq_along(pval_reduced_list)) {
      # Get significant features for this iteration (both reduced and removed sets)
      if (!use_ABC_for_feature_selection || !requireNamespace("ABCanalysis", quietly = TRUE)) {
        reduced_significant_raw <- which(pval_reduced_list[[i]]$p_values < significance_threshold)
        reduced_fdr_adjusted <- p.adjust(pval_reduced_list[[i]]$p_values, method = "fdr")
        reduced_significant_fdr <- which(reduced_fdr_adjusted < significance_threshold)

        removed_significant_raw <- which(pval_removed_list[[i]]$p_values < significance_threshold)
        removed_fdr_adjusted <- p.adjust(pval_removed_list[[i]]$p_values, method = "fdr")
        removed_significant_fdr <- which(removed_fdr_adjusted < significance_threshold)
      } else {
        reduced_significant_raw <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(pval_reduced_list[[i]]$p_values, 1e-100)))$Aind,
          paste("ABCanalysis for reduced raw p-values (iteration", i, ")"),
          which(pval_reduced_list[[i]]$p_values < significance_threshold)
        )

        reduced_fdr_adjusted <- p.adjust(pval_reduced_list[[i]]$p_values, method = "fdr")
        reduced_significant_fdr <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(reduced_fdr_adjusted, 1e-100)))$Aind,
          paste("ABCanalysis for reduced FDR-corrected p-values (iteration", i, ")"),
          which(reduced_fdr_adjusted < significance_threshold)
        )

        removed_significant_raw <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(pval_removed_list[[i]]$p_values, 1e-100)))$Aind,
          paste("ABCanalysis for removed raw p-values (iteration", i, ")"),
          which(pval_removed_list[[i]]$p_values < significance_threshold)
        )

        removed_fdr_adjusted <- p.adjust(pval_removed_list[[i]]$p_values, method = "fdr")
        removed_significant_fdr <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(removed_fdr_adjusted, 1e-100)))$Aind,
          paste("ABCanalysis for removed FDR-corrected p-values (iteration", i, ")"),
          which(removed_fdr_adjusted < significance_threshold)
        )
      }

      # Calculate Jaccard indices for all three comparisons
      jaccard_results <- rbind(jaccard_results, data.frame(
        Iteration = i,
        Comparison = "Reduced vs Original",
        Jaccard_Raw = calculate_jaccard_index(original_significant_raw, reduced_significant_raw),
        Jaccard_FDR = calculate_jaccard_index(original_significant_fdr, reduced_significant_fdr),
        stringsAsFactors = FALSE
      ))

      jaccard_results <- rbind(jaccard_results, data.frame(
        Iteration = i,
        Comparison = "Removed vs Original",
        Jaccard_Raw = calculate_jaccard_index(original_significant_raw, removed_significant_raw),
        Jaccard_FDR = calculate_jaccard_index(original_significant_fdr, removed_significant_fdr),
        stringsAsFactors = FALSE
      ))

      jaccard_results <- rbind(jaccard_results, data.frame(
        Iteration = i,
        Comparison = "Reduced vs Removed",
        Jaccard_Raw = calculate_jaccard_index(reduced_significant_raw, removed_significant_raw),
        Jaccard_FDR = calculate_jaccard_index(reduced_significant_fdr, removed_significant_fdr),
        stringsAsFactors = FALSE
      ))
    }

    return(list(
      jaccard_results = jaccard_results,
      original_significant_raw = original_significant_raw,
      original_significant_fdr = original_significant_fdr,
      fdr_adjusted_pvals = fdr_adjusted_pvals
    ))
  }

  # Create Jaccard index panel for integration with comparison plot
  create_jaccard_panel <- function(jaccard_analysis) {
    # Prepare data for plotting - reshape from wide to long format
    jaccard_long <- reshape2::melt(jaccard_analysis$jaccard_results,
                                   id.vars = c("Iteration", "Comparison"),
                                   variable.name = "Correction_Type",
                                   value.name = "Jaccard_Index")

    # Clean up labels for better readability
    jaccard_long$Correction_Type <- ifelse(jaccard_long$Correction_Type == "Jaccard_Raw",
                                           "Raw p-values", "FDR corrected")

    # Ensure proper factor ordering for consistent colors across plots
    jaccard_long$Comparison <- factor(jaccard_long$Comparison,
                                      levels = c("Reduced vs Original", "Removed vs Original", "Reduced vs Removed"))

    # Create the panel with consistent styling
    plot_theme <- theme_light() +
      theme(
        legend.position = "none",
        strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 11)
      )

    # Create Jaccard index visualization with reference lines
    ggplot(jaccard_long, aes(x = Comparison, y = Jaccard_Index, color = Comparison)) +
      geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk",
                   outlier.shape = NA, linewidth = 0.3) +
      geom_point(aes(color = Comparison, shape = as.factor(Iteration)),
                 position = position_jitter(width = 0.1), alpha = 0.7, size = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3", linewidth = 0.5) +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 0.3) +
      geom_hline(yintercept = 0, linetype = "solid", color = "red", linewidth = 0.3, alpha = 0.5) +
      facet_wrap(~Correction_Type, scales = "free_x", ncol = 1) +
      scale_color_colorblind() +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
      plot_theme +
      labs(
        x = "Comparison",
        y = "Jaccard index",
        title = "Feature selection overlap"
      ) +
      guides(shape = "none")
  }

  # Execute Jaccard analysis workflow
  cat("Calculating Jaccard indices for feature selection overlap...\n")
  jaccard_analysis <- calculate_pvals_with_jaccard_threefold(pval_results, pval_reduced_list, pval_removed_list,
                                                             significance_threshold)

  # Create the Jaccard index panel
  p_jaccard <- create_jaccard_panel(jaccard_analysis)

  # Print summary statistics for Jaccard indices
  cat("\nJaccard index summary:\n")
  cat("=================================\n")

  # Calculate summary statistics for each comparison type
  summary_stats <- aggregate(cbind(Jaccard_Raw, Jaccard_FDR) ~ Comparison,
                             data = jaccard_analysis$jaccard_results,
                             FUN = function(x) c(Mean = round(mean(x, na.rm = TRUE), 3),
                                                 Median = round(median(x, na.rm = TRUE), 3),
                                                 SD = round(sd(x, na.rm = TRUE), 3),
                                                 Min = round(min(x, na.rm = TRUE), 3),
                                                 Max = round(max(x, na.rm = TRUE), 3)))

  print(summary_stats)

  # Additional insights about feature selection
  cat("\nFeature selection insights:\n")
  cat("===========================\n")
  cat(paste("Original significant features (raw p < ", significance_threshold, "): ",
            length(jaccard_analysis$original_significant_raw), " out of ", nrow(pval_results), "\n", sep = ""))
  cat(paste("Original significant features (FDR corrected): ",
            length(jaccard_analysis$original_significant_fdr), " out of ", nrow(pval_results), "\n", sep = ""))

  # Interpretation guidelines for users
  cat("\nInterpretation guidelines:\n")
  cat("=========================\n")
  cat("• Reduced vs Original: How well does downsampling preserve significant features?\n")
  cat("• Removed vs Original: Do removed samples contain different significant features?\n")
  cat("• Reduced vs Removed: How complementary are the feature sets in the two partitions?\n")
  cat("• High Jaccard (>0.7): Strong overlap\n")
  cat("• Medium Jaccard (0.3-0.7): Moderate overlap\n")
  cat("• Low Jaccard (<0.3): Weak overlap\n")

  
  # ---------------------------------------------------------------------------
  # SIMPLE MATCHING COEFFICIENT ANALYSIS (Feature Selection Overlap)
  # ---------------------------------------------------------------------------
  
  # Simple Matching Coefficient function in R
  # vec1, vec2: vectors of equal length, values can be "yes"/"no" or 1/0
  
  simple_matching_coefficient <- function(vec1, vec2) {
    # Check inputs are the same length
    if(length(vec1) != length(vec2)) stop("Vectors must be of equal length")
    
    # Convert to binary (1 for "yes", 0 for "no") if needed
    # Skip this conversion if already using 1 and 0
    if(is.character(vec1)) vec1 <- as.integer(tolower(vec1) == "yes")
    if(is.character(vec2)) vec2 <- as.integer(tolower(vec2) == "yes")
    
    # Calculate number of matches
    matches <- sum(vec1 == vec2)
    
    # Calculate SMC as proportion of matches
    smc <- matches / length(vec1)
    return(smc)
  }
  
  
  # P-value calculation with simple matching coefficient analysis
  # Compares feature selection overlap between: Original-Reduced, Original-Removed, Reduced-Removed
  calculate_pvals_with_simple_matching_coefficient_threefold <- function(pval_results, pval_reduced_list, pval_removed_list,
                                                                         significance_threshold = 0.05) {
    
    # Apply FDR correction for multiple testing
    fdr_adjusted_pvals <- p.adjust(pval_results$p_values, method = "fdr")
    
    cat("Calculating simple matching coefficient indices for feature selection overlap...\n")
    
    # Create vector of zeros for simple matching coefficient
    simple_matching_coefficient_zeros <- rep(0, length(pval_results$Variable_No))
    
    # Identify significant features for original data using both raw and FDR-corrected p-values
    if (!use_ABC_for_feature_selection || !requireNamespace("ABCanalysis", quietly = TRUE)) {
      original_significant_raw <- which(pval_results$p_values < significance_threshold)
      original_significant_raw_01 <- simple_matching_coefficient_zeros
      original_significant_raw_01[original_significant_raw] <- 1
      original_significant_fdr <- which(fdr_adjusted_pvals < significance_threshold)
      original_significant_fdr_01 <- simple_matching_coefficient_zeros
      original_significant_fdr_01[original_significant_raw] <- 1
    } else {
      # Use ABCanalysis for feature selection based on -log10 transformed p-values
      original_significant_raw <- safe_test_execution(
        function() ABCanalysis::ABCanalysis(-log10(pmax(pval_results$p_values, 1e-100)))$Aind,
        "ABCanalysis for raw p-values",
        which(pval_results$p_values < significance_threshold)
      )
      original_significant_raw_01 <- simple_matching_coefficient_zeros
      original_significant_raw_01[original_significant_raw] <- 1
      
      original_significant_fdr <- safe_test_execution(
        function() ABCanalysis::ABCanalysis(-log10(pmax(fdr_adjusted_pvals, 1e-100)))$Aind,
        "ABCanalysis for FDR-corrected p-values",
        which(fdr_adjusted_pvals < significance_threshold)
      )
      original_significant_fdr_01 <- simple_matching_coefficient_zeros
      original_significant_fdr_01[original_significant_raw] <- 1
    }
    
    # Initialize results data frame for all three comparison types
    simple_matching_coefficient_results <- data.frame(
      Iteration = integer(),
      Comparison = character(),
      simple_matching_coefficient_Raw = numeric(),
      simple_matching_coefficient_FDR = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Process each downsampling iteration
    for (i in seq_along(pval_reduced_list)) {
      # Get significant features for this iteration (both reduced and removed sets)
      if (!use_ABC_for_feature_selection || !requireNamespace("ABCanalysis", quietly = TRUE)) {
        reduced_significant_raw <- which(pval_reduced_list[[i]]$p_values < significance_threshold)
        reduced_significant_raw_01 <- simple_matching_coefficient_zeros
        reduced_significant_raw_01[reduced_significant_raw] <- 1
        reduced_fdr_adjusted <- p.adjust(pval_reduced_list[[i]]$p_values, method = "fdr")
        reduced_significant_fdr <- which(reduced_fdr_adjusted < significance_threshold)
        reduced_significant_fdr_01 <- simple_matching_coefficient_zeros
        reduced_significant_fdr_01[reduced_significant_fdr] <- 1
        
        removed_significant_raw <- which(pval_removed_list[[i]]$p_values < significance_threshold)
        removed_significant_raw_01 <- simple_matching_coefficient_zeros
        removed_significant_raw_01[removed_significant_raw] <- 1
        removed_fdr_adjusted <- p.adjust(pval_removed_list[[i]]$p_values, method = "fdr")
        removed_significant_fdr <- which(removed_fdr_adjusted < significance_threshold)
        removed_significant_fdr_01 <- simple_matching_coefficient_zeros
        removed_significant_fdr_01[removed_significant_fdr] <- 1
      } else {
        reduced_significant_raw <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(pval_reduced_list[[i]]$p_values, 1e-100)))$Aind,
          paste("ABCanalysis for reduced raw p-values (iteration", i, ")"),
          which(pval_reduced_list[[i]]$p_values < significance_threshold)
        )
        reduced_significant_raw_01 <- simple_matching_coefficient_zeros
        reduced_significant_raw_01[reduced_significant_raw] <- 1
        
        reduced_fdr_adjusted <- p.adjust(pval_reduced_list[[i]]$p_values, method = "fdr")
        reduced_significant_fdr <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(reduced_fdr_adjusted, 1e-100)))$Aind,
          paste("ABCanalysis for reduced FDR-corrected p-values (iteration", i, ")"),
          which(reduced_fdr_adjusted < significance_threshold)
        )
        reduced_significant_fdr_01 <- simple_matching_coefficient_zeros
        reduced_significant_fdr_01[reduced_significant_fdr] <- 1
        
        removed_significant_raw <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(pval_removed_list[[i]]$p_values, 1e-100)))$Aind,
          paste("ABCanalysis for removed raw p-values (iteration", i, ")"),
          which(pval_removed_list[[i]]$p_values < significance_threshold)
        )
        removed_significant_raw_01 <- simple_matching_coefficient_zeros
        removed_significant_raw_01[removed_significant_raw] <- 1
        
        removed_fdr_adjusted <- p.adjust(pval_removed_list[[i]]$p_values, method = "fdr")
        removed_significant_fdr <- safe_test_execution(
          function() ABCanalysis::ABCanalysis(-log10(pmax(removed_fdr_adjusted, 1e-100)))$Aind,
          paste("ABCanalysis for removed FDR-corrected p-values (iteration", i, ")"),
          which(removed_fdr_adjusted < significance_threshold)
        )
        removed_significant_fdr_01 <- simple_matching_coefficient_zeros
        removed_significant_fdr_01[removed_significant_fdr] <- 1
      }
      
      # Calculate simple matching coefficient indices for all three comparisons
      simple_matching_coefficient_results <- rbind(simple_matching_coefficient_results, data.frame(
        Iteration = i,
        Comparison = "Reduced vs Original",
        simple_matching_coefficient_Raw = simple_matching_coefficient(original_significant_raw_01, reduced_significant_raw_01),
        simple_matching_coefficient_FDR = simple_matching_coefficient(original_significant_fdr_01, reduced_significant_fdr_01),
        stringsAsFactors = FALSE
      ))
      
      simple_matching_coefficient_results <- rbind(simple_matching_coefficient_results, data.frame(
        Iteration = i,
        Comparison = "Removed vs Original",
        simple_matching_coefficient_Raw = simple_matching_coefficient(original_significant_raw_01, removed_significant_raw_01),
        simple_matching_coefficient_FDR = simple_matching_coefficient(original_significant_fdr_01, removed_significant_fdr_01),
        stringsAsFactors = FALSE
      ))
      
      simple_matching_coefficient_results <- rbind(simple_matching_coefficient_results, data.frame(
        Iteration = i,
        Comparison = "Reduced vs Removed",
        simple_matching_coefficient_Raw = simple_matching_coefficient(reduced_significant_raw_01, removed_significant_raw_01),
        simple_matching_coefficient_FDR = simple_matching_coefficient(reduced_significant_fdr_01, removed_significant_fdr_01),
        stringsAsFactors = FALSE
      ))
    }
    
    return(list(
      simple_matching_coefficient_results = simple_matching_coefficient_results,
      original_significant_raw_01 = original_significant_raw_01,
      original_significant_fdr_01 = original_significant_fdr_01,
      fdr_adjusted_pvals = fdr_adjusted_pvals
    ))
  }
  
  # Create simple matching coefficient panel for integration with comparison plot
  create_simple_matching_coefficient_panel <- function(simple_matching_coefficient_analysis) {
    # Prepare data for plotting - reshape from wide to long format
    simple_matching_coefficient_long <- reshape2::melt(simple_matching_coefficient_analysis$simple_matching_coefficient_results,
                                                       id.vars = c("Iteration", "Comparison"),
                                                       variable.name = "Correction_Type",
                                                       value.name = "simple_matching_coefficient")
    
    # Clean up labels for better readability
    simple_matching_coefficient_long$Correction_Type <- ifelse(simple_matching_coefficient_long$Correction_Type == "simple_matching_coefficient_Raw",
                                                               "Raw p-values", "FDR corrected")
    
    # Ensure proper factor ordering for consistent colors across plots
    simple_matching_coefficient_long$Comparison <- factor(simple_matching_coefficient_long$Comparison,
                                                          levels = c("Reduced vs Original", "Removed vs Original", "Reduced vs Removed"))
    
    # Create the panel with consistent styling
    plot_theme <- theme_light() +
      theme(
        legend.position = "none",
        strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 11)
      )
    
    # Create simple matching coefficient visualization with reference lines
    ggplot(simple_matching_coefficient_long, aes(x = Comparison, y = simple_matching_coefficient, color = Comparison)) +
      geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk",
                   outlier.shape = NA, linewidth = 0.3) +
      geom_point(aes(color = Comparison, shape = as.factor(Iteration)),
                 position = position_jitter(width = 0.1), alpha = 0.7, size = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3", linewidth = 0.5) +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 0.3) +
      geom_hline(yintercept = 0, linetype = "solid", color = "red", linewidth = 0.3, alpha = 0.5) +
      facet_wrap(~Correction_Type, scales = "free_x", ncol = 1) +
      scale_color_colorblind() +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
      plot_theme +
      labs(
        x = "Comparison",
        y = "Simple matching coefficient",
        title = "Feature selection overlap"
      ) +
      guides(shape = "none")
  }
  
  # Execute simple matching coefficient analysis workflow
  cat("Calculating simple matching coefficient indices for feature selection overlap...\n")
  simple_matching_coefficient_analysis <- calculate_pvals_with_simple_matching_coefficient_threefold(pval_results, pval_reduced_list, pval_removed_list,
                                                                                                     significance_threshold)
  
  # Create the simple matching coefficient panel
  p_simple_matching_coefficient <- create_simple_matching_coefficient_panel(simple_matching_coefficient_analysis)
  
  # Print summary statistics for Jaccard indices
  cat("\nSimple matching coefficient Summary:\n")
  cat("=================================\n")
  
  # Calculate summary statistics for each comparison type
  summary_stats <- aggregate(cbind(simple_matching_coefficient_Raw, simple_matching_coefficient_FDR) ~ Comparison,
                             data = simple_matching_coefficient_analysis$simple_matching_coefficient_results,
                             FUN = function(x) c(Mean = round(mean(x, na.rm = TRUE), 3),
                                                 Median = round(median(x, na.rm = TRUE), 3),
                                                 SD = round(sd(x, na.rm = TRUE), 3),
                                                 Min = round(min(x, na.rm = TRUE), 3),
                                                 Max = round(max(x, na.rm = TRUE), 3)))
  
  print(summary_stats)
  
  # Additional insights about feature selection
  cat("\nFeature Selection Insights:\n")
  cat("===========================\n")
  cat(paste("Original significant features (raw p < ", significance_threshold, "): ",
            length(simple_matching_coefficient_analysis$original_significant_raw), " out of ", nrow(pval_results), "\n", sep = ""))
  cat(paste("Original significant features (FDR corrected): ",
            length(simple_matching_coefficient_analysis$original_significant_fdr), " out of ", nrow(pval_results), "\n", sep = ""))
  
  # Interpretation guidelines for users
  cat("\nInterpretation Guidelines:\n")
  cat("=========================\n")
  cat("• Reduced vs Original: How well does downsampling preserve significant features?\n")
  cat("• Removed vs Original: Do removed samples contain different significant features?\n")
  cat("• Reduced vs Removed: How complementary are the feature sets in the two partitions?\n")
  cat("• High simple matching coefficient (>0.7): Strong overlap\n")
  cat("• Medium simple matching coefficient (0.3-0.7): Moderate overlap\n")
  cat("• Low simple matching coefficient (<0.3): Weak overlap\n")
  
  # ---------------------------------------------------------------------------
  # CORRELATION ANALYSIS (Kendall's Tau)
  # ---------------------------------------------------------------------------

  cat("Performing correlation analysis...\n")

  # Validation and error handling
  if (length(pval_reduced_list) == 0 || length(pval_removed_list) == 0) {
    stop("Empty p-value lists provided")
  }

  # Extract p-values matrices for correlation analysis
  pval_reduced_matrix <- do.call(cbind, lapply(pval_reduced_list, `[[`, "p_values"))
  pval_removed_matrix <- do.call(cbind, lapply(pval_removed_list, `[[`, "p_values"))

  # Validate matrices
  if (!is.numeric(pval_reduced_matrix) || !is.numeric(pval_removed_matrix)) {
    stop("Non-numeric p-values detected")
  }
  if (any(is.na(pval_reduced_matrix)) || any(is.na(pval_removed_matrix))) {
    warning("NA values detected in p-value matrices")
  }

  # Original p-values
  pval_orig <- pval_results$p_values
  if (is.null(pval_orig) || !is.numeric(pval_orig)) {
    stop("Invalid original p-values")
  }

  # Filter features for the originally relevant features
  if (use_ABC_for_Pvalues && requireNamespace("ABCanalysis", quietly = TRUE)) {
    feature_indices_orig_C <- safe_test_execution(
      function() ABCanalysis::ABCanalysis(-log10(pmax(pval_orig, 1e-100)))$Cind,
      "ABCanalysis for feature filtering",
      integer(0)
    )

    if (length(feature_indices_orig_C) == 0) {
      cat("No features selected by ABCanalysis - using all features\n")
    } else {
      # Validate indices
      valid_indices <- feature_indices_orig_C[feature_indices_orig_C <= length(pval_orig) &
                                                feature_indices_orig_C >= 1]

      if (length(valid_indices) < length(feature_indices_orig_C)) {
        warning("Some invalid feature indices detected, filtering them out")
        feature_indices_orig_C <- valid_indices
      }

      # Check remaining dimensions after filtering
      remaining_features <- length(pval_orig) - length(feature_indices_orig_C)
      if (remaining_features <= 1) {
        cat("Too few features remaining after ABC filtering - using all features\n")
      } else {
        # Apply filtering
        pval_orig <- pval_orig[-feature_indices_orig_C]
        pval_removed_matrix <- pval_removed_matrix[-feature_indices_orig_C,, drop = FALSE]
        pval_reduced_matrix <- pval_reduced_matrix[-feature_indices_orig_C,, drop = FALSE]
        cat(paste("ABC filtering removed", length(feature_indices_orig_C), "features\n"))
      }
    }
  }

  # Ensure pval_orig is a vector for correlation analysis
  pval_orig <- as.vector(pval_orig)

  # Compute Kendall's tau correlations with error handling
  tau_reduced <- apply(pval_reduced_matrix, 2, function(x) {
    safe_test_execution(
      function() cor(x, pval_orig, method = "kendall", use = "complete.obs"),
      "Kendall's tau correlation for reduced vs original"
    )
  })

  tau_removed <- apply(pval_removed_matrix, 2, function(x) {
    safe_test_execution(
      function() cor(x, pval_orig, method = "kendall", use = "complete.obs"),
      "Kendall's tau correlation for removed vs original"
    )
  })

  tau_reduced_vs_removed <- mapply(function(x, y) {
    safe_test_execution(
      function() cor(x, y, method = "kendall", use = "complete.obs"),
      "Kendall's tau correlation for reduced vs removed"
    )
  }, as.data.frame(pval_reduced_matrix), as.data.frame(pval_removed_matrix))

  # Validate correlation results
  if (all(is.na(tau_reduced)) || all(is.na(tau_removed)) || all(is.na(tau_reduced_vs_removed))) {
    warning("All correlations are NA - check data quality")
  }

  # Prepare data for plotting correlations
  correlation_data <- rbind(
    data.frame(Iteration = seq_along(tau_reduced), Tau = tau_reduced, Comparison = "Reduced vs Original"),
    data.frame(Iteration = seq_along(tau_removed), Tau = tau_removed, Comparison = "Removed vs Original"),
    data.frame(Iteration = seq_along(tau_reduced_vs_removed), Tau = tau_reduced_vs_removed, Comparison = "Reduced vs Removed")
  )

  # Calculate summary statistics for annotations with NA handling
  calc_safe_stat <- function(x, fun) {
    if (all(is.na(x))) return(NA)

    # Special handling for mode function
    if (identical(fun, mode)) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return(NA)

      # For continuous data, we'll use density-based mode estimation
      # Find the value with highest density
      density_result <- density(x, na.rm = TRUE)
      mode_value <- density_result$x[which.max(density_result$y)]
      return(mode_value)
    }

    # For other functions, use na.rm = TRUE
    fun(x, na.rm = TRUE)
  }

  median_tau_reduced <- calc_safe_stat(tau_reduced, median)
  median_tau_removed <- calc_safe_stat(tau_removed, median)
  median_tau_reduced_vs_removed <- calc_safe_stat(tau_reduced_vs_removed, median)

  variance_tau_reduced <- calc_safe_stat(tau_reduced, var)
  variance_tau_removed <- calc_safe_stat(tau_removed, var)
  variance_tau_reduced_vs_removed <- calc_safe_stat(tau_reduced_vs_removed, var)

  min_tau_reduced <- calc_safe_stat(tau_reduced, min)
  min_tau_removed <- calc_safe_stat(tau_removed, min)
  min_tau_reduced_vs_removed <- calc_safe_stat(tau_reduced_vs_removed, min)

  mode_tau_reduced <- calc_safe_stat(tau_reduced, mode)
  mode_tau_removed <- calc_safe_stat(tau_removed, mode)
  mode_tau_reduced_vs_removed <- calc_safe_stat(tau_reduced_vs_removed, mode)

  # Create annotation data frame for plot statistics - 2x2 matrix arrangement
  annotation_df <- data.frame(
    Comparison = rep(c("Reduced vs Original", "Removed vs Original", "Reduced vs Removed"), 4),
    x = rep(c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2), 2), # Side by side positioning for each group
    y = c(rep(1.25, 6), rep(1.15, 6)), # Two rows
    color_group = rep(c("Reduced vs Original", "Reduced vs Original",
                        "Reduced vs Removed", "Reduced vs Removed",
                        "Removed vs Original", "Removed vs Original"), 2),
    label = c(
  # Top row - Med and Var values side by side
      paste("Med:", round(median_tau_reduced, 3)),
      paste("Var:", round(variance_tau_reduced, 3)),
      paste("Med:", round(median_tau_reduced_vs_removed, 3)),
      paste("Var:", round(variance_tau_reduced_vs_removed, 3)),
      paste("Med:", round(median_tau_removed, 3)),
      paste("Var:", round(variance_tau_removed, 3)),
  # Bottom row - Min and Mode values side by side
      paste("Min:", round(min_tau_reduced, 3)),
      paste("Mode:", round(mode_tau_reduced, 3)),
      paste("Min:", round(min_tau_reduced_vs_removed, 3)),
      paste("Mode:", round(mode_tau_reduced_vs_removed, 3)),
      paste("Min:", round(min_tau_removed, 3)),
      paste("Mode:", round(mode_tau_removed, 3))
    )
  )

  # Create Kendall's tau correlation plot with box plots and individual points
  p_kendall_correlations_box <- ggplot(correlation_data, aes(x = factor(Comparison), y = Tau)) +
    geom_boxplot(aes(color = Comparison), position = "dodge", alpha = 0.5, fill = "cornsilk") +
    geom_point(aes(color = Comparison, shape = as.factor(Iteration)), position = position_jitter(width = 0.05), alpha = 1) +
    ggthemes::scale_color_colorblind() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3") +
    geom_text(
      data = annotation_df,
      aes(x = x, y = y, label = label, color = color_group),
      inherit.aes = FALSE,
      hjust = 1, vjust = 0.5, size = 2.5, angle = 90, fontface = "plain"
    ) +
    theme_light() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      x = "Comparison",
      y = "Kendall's Tau",
      title = "Kendall's Tau",
      subtitle = "p (group differences)"
    ) +
    guides(shape = "none")

  # ---------------------------------------------------------------------------
  # FINAL VISUALIZATION AND OUTPUT
  # ---------------------------------------------------------------------------

  # Combine and display all plots in a comprehensive layout
  combined_plot <- cowplot::plot_grid(comparison_plot, 
                                      p_kendall_correlations_box, 
                                      p_jaccard,
                                      p_simple_matching_coefficient,
                                      labels = "AUTO", rel_widths = c(4, 1, 1, 1), nrow = 1)
  print(combined_plot)

  # Create safe filename for saving
  safe_filename <- paste0("combined_plot_",
                          format(nTrials, scientific = FALSE), "trials_per_iteration_",
                          downsampling_size, "sampled", "_CheckRemoved", CheckRemoved,
                          "_CheckThreefold", CheckThreefold,
                          "_OptimizeBetween", OptimizeBetween,
                          "_NonNoiseSelection", NonNoiseSelection, ".svg")

  # Save combined plot to file with error handling
  tryCatch({
    ggsave(safe_filename, combined_plot, width = 12, height = 12)
    cat(paste("Plot saved as:", safe_filename, "\n"))
  }, error = function(e) {
    warning(paste("Failed to save plot:", e$message), call. = FALSE)
  })

  # Print summary statistics to console for quick reference
  print(paste("Reduced vs Original - Med:", round(median_tau_reduced, 3),
              "Var:", round(variance_tau_reduced, 3),
              "Min:", round(min_tau_reduced, 3),
              "Mode:", round(mode_tau_reduced, 3)))
  print(paste("Removed vs Original - Med:", round(median_tau_removed, 3),
              "Var:", round(variance_tau_removed, 3),
              "Min:", round(min_tau_removed, 3),
              "Mode:", round(mode_tau_removed, 3)))
  print(paste("Reduced vs Removed - Med:", round(median_tau_reduced_vs_removed, 3),
              "Var:", round(variance_tau_reduced_vs_removed, 3),
              "Min:", round(min_tau_reduced_vs_removed, 3),
              "Mode:", round(mode_tau_reduced_vs_removed, 3)))

  cat("\nAnalysis complete for current parameter set.\n\n")

  # Return comprehensive results for further analysis
  return(list(
  # Core results
    p_values_orig = pval_results,
    pval_reduced_list = pval_reduced_list,
    pval_removed_list = pval_removed_list,
    downsampled_data = downsampled_results,
    correlations = correlation_data,
    jaccard_indices = jaccard_analysis$jaccard_results,
    simple_matching_coefficients = simple_matching_coefficient_analysis$simple_matching_coefficient_results,

  # Plots
    combined_plot = combined_plot,
    comparison_plot = comparison_plot,
    p_kendall_correlations_box = p_kendall_correlations_box,
    p_jaccard = p_jaccard,
    p_simple_matching_coefficient = p_simple_matching_coefficient,
  

  # Analysis results
    jaccard_analysis = jaccard_analysis,
    simple_matching_coefficient_analysis = simple_matching_coefficient_analysis,
    redundancy_found = redundancy_found,

  # Summary statistics
    tau_summary = list(
      reduced_vs_original = list(median = median_tau_reduced, variance = variance_tau_reduced, min = min_tau_reduced),
      removed_vs_original = list(median = median_tau_removed, variance = variance_tau_removed, min = min_tau_removed),
      reduced_vs_removed = list(median = median_tau_reduced_vs_removed, variance = variance_tau_reduced_vs_removed, min = min_tau_reduced_vs_removed)
    ),

  # Configuration
    config = list(
      CheckRemoved = CheckRemoved,
      CheckThreefold = CheckThreefold,
      OptimizeBetween = OptimizeBetween,
      NonNoiseSelection = NonNoiseSelection,
      nSamples = length(downsampled_results),
      nTrials = nTrials,
      TestStat = TestStat,
      significance_threshold = significance_threshold,
      use_ABC_for_Pvalues = use_ABC_for_Pvalues,
      use_ABC_for_feature_selection = use_ABC_for_feature_selection,
      use_regression_for_p = use_regression_for_p
    )
  ))
}


# ===============================================================================
# POST-PROCESSING AND COMPARATIVE ANALYSIS
# ===============================================================================

#' Extract and process correlation data from experiment results
#' @param results_single Single trial results (optional)
#' @param results_multiple Multiple trial results (optional) 
#' @param default_param_name Default parameter combination name for single trial
#' @return Combined correlation data frame with parameter information
extract_correlation_data <- function(results_single = NULL, results_multiple = NULL,
                                     default_param_name = "OpF_NoF") {
  
  combined_correlations <- data.frame()
  
  # Process single trial results if available
  if (!is.null(results_single) && !is.null(results_single$correlations)) {
    single_correlations <- results_single$correlations
    single_correlations$ParameterCombination <- default_param_name
    single_correlations$SingleMultiple <- "Single"
    combined_correlations <- rbind(combined_correlations, single_correlations)
    
    if (exists("cat")) {
      cat("Extracted correlation data from single trial experiment\n")
    }
  }
  
  # Process multiple trial results if available
  if (!is.null(results_multiple) && is.list(results_multiple)) {
    multiple_correlations <- do.call(rbind, lapply(names(results_multiple), function(param_name) {
      if (!is.null(results_multiple[[param_name]]$correlations)) {
        correlations_df <- results_multiple[[param_name]]$correlations
        correlations_df$ParameterCombination <- param_name
        return(correlations_df)
      }
      return(NULL)
    }))
    
    if (!is.null(multiple_correlations)) {
      rownames(multiple_correlations) <- NULL
      multiple_correlations$SingleMultiple <- "Multiple"
      combined_correlations <- rbind(combined_correlations, multiple_correlations)
      
      if (exists("cat")) {
        cat(sprintf("Extracted correlation data from %d parameter combinations\n",
                    length(unique(multiple_correlations$ParameterCombination))))
      }
    }
  }
  
  # Parse parameter combinations into separate columns if data exists
  if (nrow(combined_correlations) > 0) {
    combined_correlations <- combined_correlations %>%
      separate(ParameterCombination,
               into = c("OptimizeBetween", "NonNoiseSelection"),
               sep = "_",
               remove = FALSE) %>%
      mutate(
        OptimizeBetween = gsub("OB", "", OptimizeBetween),
        NonNoiseSelection = gsub("NNS", "", NonNoiseSelection)
      )
  }
  
  return(combined_correlations)
}

#' Extract and process Jaccard index data from experiment results
#' @param results_single Single trial results (optional)
#' @param results_multiple Multiple trial results (optional)
#' @param default_param_name Default parameter combination name for single trial
#' @return Combined Jaccard data frame with parameter information
extract_jaccard_data <- function(results_single = NULL, results_multiple = NULL,
                                 default_param_name = "OpF_NoF") {
  
  combined_jaccard <- data.frame()
  
  # Process single trial results if available
  if (!is.null(results_single) && !is.null(results_single$jaccard_indices)) {
    single_jaccard <- results_single$jaccard_indices
    single_jaccard$ParameterCombination <- default_param_name
    single_jaccard$SingleMultiple <- "Single"
    combined_jaccard <- rbind(combined_jaccard, single_jaccard)
    
    if (exists("cat")) {
      cat("Extracted Jaccard data from single trial experiment\n")
    }
  }
  
  # Process multiple trial results if available
  if (!is.null(results_multiple) && is.list(results_multiple)) {
    multiple_jaccard <- do.call(rbind, lapply(names(results_multiple), function(param_name) {
      if (!is.null(results_multiple[[param_name]]$jaccard_indices)) {
        jaccard_df <- results_multiple[[param_name]]$jaccard_indices
        jaccard_df$ParameterCombination <- param_name
        return(jaccard_df)
      }
      return(NULL)
    }))
    
    if (!is.null(multiple_jaccard)) {
      rownames(multiple_jaccard) <- NULL
      multiple_jaccard$SingleMultiple <- "Multiple"
      combined_jaccard <- rbind(combined_jaccard, multiple_jaccard)
      
      if (exists("cat")) {
        cat(sprintf("Extracted Jaccard data from %d parameter combinations\n",
                    length(unique(multiple_jaccard$ParameterCombination))))
      }
    }
  }
  
  # Parse parameter combinations into separate columns if data exists
  if (nrow(combined_jaccard) > 0) {
    combined_jaccard <- combined_jaccard %>%
      separate(ParameterCombination,
               into = c("OptimizeBetween", "NonNoiseSelection"),
               sep = "_",
               remove = FALSE) %>%
      mutate(
        OptimizeBetween = gsub("OB", "", OptimizeBetween),
        NonNoiseSelection = gsub("NNS", "", NonNoiseSelection)
      )
  }
  
  return(combined_jaccard)
}

#' Extract and process simple matching coefficient index data from experiment results
#' @param results_single Single trial results (optional)
#' @param results_multiple Multiple trial results (optional)
#' @param default_param_name Default parameter combination name for single trial
#' @return Combined simple matching coefficient data frame with parameter information
extract_simple_matching_coefficient_data <- function(results_single = NULL, results_multiple = NULL,
                                                     default_param_name = "OpF_NoF") {
  
  combined_simple_matching_coefficient <- data.frame()
  
  # Process single trial results if available
  if (!is.null(results_single) && !is.null(results_single$simple_matching_coefficients)) {
    single_simple_matching_coefficient <- results_single$simple_matching_coefficients
    single_simple_matching_coefficient$ParameterCombination <- default_param_name
    single_simple_matching_coefficient$SingleMultiple <- "Single"
    combined_simple_matching_coefficient <- rbind(combined_simple_matching_coefficient, single_simple_matching_coefficient)
    
    if (exists("cat")) {
      cat("Extracted simple matching coefficient data from single trial experiment\n")
    }
  }
  
  # Process multiple trial results if available
  if (!is.null(results_multiple) && is.list(results_multiple)) {
    multiple_simple_matching_coefficient <- do.call(rbind, lapply(names(results_multiple), function(param_name) {
      if (!is.null(results_multiple[[param_name]]$simple_matching_coefficients)) {
        simple_matching_coefficient_df <- results_multiple[[param_name]]$simple_matching_coefficients
        simple_matching_coefficient_df$ParameterCombination <- param_name
        return(simple_matching_coefficient_df)
      }
      return(NULL)
    }))
    
    if (!is.null(multiple_simple_matching_coefficient)) {
      rownames(multiple_simple_matching_coefficient) <- NULL
      multiple_simple_matching_coefficient$SingleMultiple <- "Multiple"
      combined_simple_matching_coefficient <- rbind(combined_simple_matching_coefficient, multiple_simple_matching_coefficient)
      
      if (exists("cat")) {
        cat(sprintf("Extracted simple matching coefficient data from %d parameter combinations\n",
                    length(unique(multiple_simple_matching_coefficient$ParameterCombination))))
      }
    }
  }
  
  # Parse parameter combinations into separate columns if data exists
  if (nrow(combined_simple_matching_coefficient) > 0) {
    combined_simple_matching_coefficient <- combined_simple_matching_coefficient %>%
      separate(ParameterCombination,
               into = c("OptimizeBetween", "NonNoiseSelection"),
               sep = "_",
               remove = FALSE) %>%
      mutate(
        OptimizeBetween = gsub("OB", "", OptimizeBetween),
        NonNoiseSelection = gsub("NNS", "", NonNoiseSelection)
      )
  }
  
  return(combined_simple_matching_coefficient)
}

#' Create correlation comparison visualization with statistical annotations
#' @param correlation_data Data frame with correlation results
#' @param plot_title Title for the plot
#' @param plot_subtitle Subtitle for the plot
#' @param facet_by_param Whether to create facets by parameter combination
#' @param add_annotations Whether to add statistical annotations (Med, Var, Min, Mode)
#' @return ggplot object
create_correlation_plot <- function(correlation_data, plot_title = "Correlation Analysis",
                                    plot_subtitle = "", facet_by_param = TRUE,
                                    add_annotations = TRUE) {
  
  if (nrow(correlation_data) == 0) {
    return(ggplot() +
             theme_void() +
             labs(title = "No correlation data available"))
  }
  
  # Calculate mode function
  calculate_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Base plot
  p <- ggplot(data = correlation_data, aes(x = Comparison, y = Tau)) +
    geom_boxplot(aes(color = Comparison), position = "dodge", alpha = 0.5,
                 fill = "cornsilk", show.legend = FALSE) +
    geom_point(aes(shape = as.factor(Iteration), color = Comparison),
               position = position_jitter(width = 0.05), alpha = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3") +
    scale_color_colorblind() +
    scale_shape_manual(values = 0:19) +
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = alpha("white", 0.5)),
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      x = "Comparison",
      y = "Kendall's Tau",
      title = plot_title,
      subtitle = plot_subtitle,
      shape = "Iteration"
    )
  
  # Add annotations if requested
  if (add_annotations) {
    # Calculate statistics for each comparison type AND parameter combination
    stats_by_comparison <- correlation_data %>%
      group_by(Comparison, ParameterCombination) %>%
      summarise(
        median_val = median(Tau, na.rm = TRUE),
        variance_val = var(Tau, na.rm = TRUE),
        min_val = min(Tau, na.rm = TRUE),
        mode_val = calculate_mode(Tau),
        .groups = "drop"
      )
    
    # Get unique comparisons in the order they appear in the plot
    unique_comparisons <- levels(factor(correlation_data$Comparison))
    if (is.null(unique_comparisons)) {
      unique_comparisons <- unique(correlation_data$Comparison)
    }
    
    # Create annotation data frame - 2x2 matrix arrangement for each comparison and parameter combination
    annotation_df <- data.frame()
    
    for (i in seq_along(unique_comparisons)) {
      comp <- unique_comparisons[i]
      
      # Get stats for this comparison across all parameter combinations
      comp_stats <- stats_by_comparison[stats_by_comparison$Comparison == comp,]
      
      if (nrow(comp_stats) > 0) {
        # Create annotations for this comparison for EACH parameter combination
        for (j in 1:nrow(comp_stats)) {
          param_combo <- comp_stats$ParameterCombination[j]
          stats <- comp_stats[j,]
          
          comp_annotations <- data.frame(
            Comparison = rep(comp, 4),
            ParameterCombination = rep(param_combo, 4), # Include parameter combination
            x = c(i - 0.1, i + 0.1, i - 0.1, i + 0.1), # Side by side positioning
            y = c(1.5, 1.5, 1.15, 1.15), # Two rows
            color_group = rep(comp, 4),
            label = c(
              paste("Med:", round(stats$median_val, 3)),
              paste("Var:", round(stats$variance_val, 3)),
              paste("Min:", round(stats$min_val, 3)),
              paste("Mode:", round(stats$mode_val, 3))
            )
          )
          annotation_df <- rbind(annotation_df, comp_annotations)
        }
      }
    }
    
    # Add annotations to plot
    if (nrow(annotation_df) > 0) {
      p <- p + geom_text(
        data = annotation_df,
        aes(x = x, y = y, label = label, color = color_group),
        inherit.aes = FALSE,
        hjust = 1, vjust = 0.5, size = 2.5, angle = 90, fontface = "plain"
      )
    }
  }
  
  # Add faceting if requested and multiple parameter combinations exist
  if (facet_by_param && length(unique(correlation_data$ParameterCombination)) > 1) {
    p <- p + facet_wrap(. ~ ParameterCombination, nrow = 1) +
      guides(shape = guide_legend(nrow = 1), color = "none")
  } else {
    p <- p + guides(color = "none", shape = "none")
  }
  
  # Fix the y-axis scaling issue
  if (add_annotations) {
    p <- p + scale_y_continuous(
      limits = c(floor(min(correlation_data$Tau, na.rm = TRUE) * 4) / 4, 1.55), # Extended upper limit for annotations
      breaks = seq(from = floor(min(correlation_data$Tau, na.rm = TRUE) * 4) / 4,
                   to = 1, by = 0.25)
    )
  } else {
    p <- p + scale_y_continuous(
      limits = c(min(correlation_data$Tau, na.rm = TRUE), 1), # Normal limits
      breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)
    )
  }
  
  return(p)
}

#' Create Jaccard index comparison visualization with statistical annotations
#' @param jaccard_data Data frame with Jaccard index results
#' @param plot_title Title for the plot
#' @param plot_subtitle Subtitle for the plot
#' @param facet_by_param Whether to create facets by parameter combination
#' @param add_annotations Whether to add statistical annotations (Med, Var, Min, Mode)
#' @return ggplot object
create_jaccard_plot <- function(jaccard_data, plot_title = "Jaccard index analysis",
                                plot_subtitle = "", facet_by_param = TRUE,
                                add_annotations = TRUE) {
  
  if (nrow(jaccard_data) == 0) {
    return(ggplot() +
             theme_void() +
             labs(title = "No Jaccard data available"))
  }
  
  # Calculate mode function
  calculate_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Reshape data for visualization
  jaccard_long <- reshape2::melt(jaccard_data,
                                 id.vars = c("Iteration", "Comparison", "ParameterCombination",
                                             "OptimizeBetween", "NonNoiseSelection"),
                                 measure.vars = c("Jaccard_Raw", "Jaccard_FDR"),
                                 variable.name = "Correction_Type",
                                 value.name = "Jaccard_Index")
  
  # Clean up labels for better readability
  jaccard_long$Correction_Type <- ifelse(jaccard_long$Correction_Type == "Jaccard_Raw",
                                         "Raw p-values", "FDR corrected")
  
  # Ensure proper factor ordering for consistent colors
  jaccard_long$Comparison <- factor(jaccard_long$Comparison,
                                    levels = c("Reduced vs Original", "Removed vs Original",
                                               "Reduced vs Removed"))
  
  # Create base plot
  p <- ggplot(jaccard_long, aes(x = Comparison, y = Jaccard_Index, color = Comparison)) +
    geom_boxplot(position = "dodge", alpha = 0.3, fill = "cornsilk",
                 outlier.shape = NA, size = 0.3) +
    geom_point(aes(shape = as.factor(Iteration)),
               position = position_jitter(width = 0.1), alpha = 0.7, size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3", linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", linewidth = 0.3, alpha = 0.5) +
    scale_color_colorblind() +
    scale_shape_manual(values = 0:19) +
    theme_light() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black", size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11),
      strip.text.y = element_text(angle = 90),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      x = "Comparison",
      y = "Jaccard index",
      title = plot_title,
      subtitle = plot_subtitle
    ) +
    guides(shape = "none", color = "none")
  
  # Add annotations if requested
  if (add_annotations) {
    # Calculate statistics for each comparison type, correction type, AND parameter combination
    stats_by_group <- jaccard_long %>%
      group_by(Comparison, Correction_Type, ParameterCombination) %>%
      summarise(
        median_val = median(Jaccard_Index, na.rm = TRUE),
        variance_val = var(Jaccard_Index, na.rm = TRUE),
        min_val = min(Jaccard_Index, na.rm = TRUE),
        mode_val = calculate_mode(Jaccard_Index),
        .groups = "drop"
      )
    
    # Get unique comparisons in the order they appear in the plot
    unique_comparisons <- levels(jaccard_long$Comparison)
    if (is.null(unique_comparisons)) {
      unique_comparisons <- unique(jaccard_long$Comparison)
    }
    
    # Create annotation data frame for each comparison, correction type, and parameter combination
    annotation_df <- data.frame()
    
    for (i in seq_along(unique_comparisons)) {
      comp <- unique_comparisons[i]
      
      # Get stats for this comparison across correction types and parameter combinations
      comp_stats <- stats_by_group[stats_by_group$Comparison == comp,]
      
      if (nrow(comp_stats) > 0) {
        # Create annotations for this comparison for EACH correction type and parameter combination
        for (j in 1:nrow(comp_stats)) {
          correction_type <- comp_stats$Correction_Type[j]
          param_combo <- comp_stats$ParameterCombination[j]
          stats <- comp_stats[j,]
          
          comp_annotations <- data.frame(
            Comparison = rep(comp, 4),
            Correction_Type = rep(correction_type, 4),
            ParameterCombination = rep(param_combo, 4), # Include parameter combination
            x = c(i - 0.1, i + 0.1, i - 0.1, i + 0.1), # Side by side positioning
            y = c(1.5, 1.5, 1.15, 1.15), # Two rows
            color_group = rep(comp, 4),
            label = c(
              paste("Med:", round(stats$median_val, 3)),
              paste("Var:", round(stats$variance_val, 3)),
              paste("Min:", round(stats$min_val, 3)),
              paste("Mode:", round(stats$mode_val, 3))
            )
          )
          annotation_df <- rbind(annotation_df, comp_annotations)
        }
      }
    }
    
    # Add annotations to plot
    if (nrow(annotation_df) > 0) {
      p <- p + geom_text(
        data = annotation_df,
        aes(x = x, y = y, label = label, color = color_group),
        inherit.aes = FALSE,
        hjust = 1, vjust = 0.5, size = 2.2, angle = 90, fontface = "plain"
      )
    }
    
    # Adjust y-axis limits to accommodate annotations
    p <- p + scale_y_continuous(limits = c(0, 1.55), breaks = seq(0, 1, 0.25))
  } else {
    # Original y-axis limits without annotations
    p <- p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))
  }
  
  # Add faceting if requested and data supports it
  if (facet_by_param) {
    if (length(unique(jaccard_long$ParameterCombination)) > 1) {
      p <- p + facet_grid(Correction_Type ~ ParameterCombination, scales = "free_x") +
        labs(caption = "Rows show p-value correction type; Columns show parameter combinations")
    } else {
      # Single parameter combination - just facet by correction type
      p <- p + facet_wrap(Correction_Type ~ ., ncol = 1, strip.position = "right")
    }
  }
  
  return(p)
}

#' Create simple matching coefficient index comparison visualization with statistical annotations
#' @param simple_matching_coefficient_data Data frame with simple matching coefficient index results
#' @param plot_title Title for the plot
#' @param plot_subtitle Subtitle for the plot
#' @param facet_by_param Whether to create facets by parameter combination
#' @param add_annotations Whether to add statistical annotations (Med, Var, Min, Mode)
#' @return ggplot object
create_simple_matching_coefficient_plot <- function(simple_matching_coefficient_data, plot_title = "Simple matching coefficient analysis",
                                                    plot_subtitle = "", facet_by_param = TRUE,
                                                    add_annotations = TRUE) {
  
  if (nrow(simple_matching_coefficient_data) == 0) {
    return(ggplot() +
             theme_void() +
             labs(title = "No simple matching coefficient data available"))
  }
  
  # Calculate mode function
  calculate_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Reshape data for visualization
  simple_matching_coefficient_long <- reshape2::melt(simple_matching_coefficient_data,
                                                     id.vars = c("Iteration", "Comparison", "ParameterCombination",
                                                                 "OptimizeBetween", "NonNoiseSelection"),
                                                     measure.vars = c("simple_matching_coefficient_Raw", "simple_matching_coefficient_FDR"),
                                                     variable.name = "Correction_Type",
                                                     value.name = "simple_matching_coefficient_Index")
  
  # Clean up labels for better readability
  simple_matching_coefficient_long$Correction_Type <- ifelse(simple_matching_coefficient_long$Correction_Type == "simple_matching_coefficient_Raw",
                                                             "Raw p-values", "FDR corrected")
  
  # Ensure proper factor ordering for consistent colors
  simple_matching_coefficient_long$Comparison <- factor(simple_matching_coefficient_long$Comparison,
                                                        levels = c("Reduced vs Original", "Removed vs Original",
                                                                   "Reduced vs Removed"))
  
  # Create base plot
  p <- ggplot(simple_matching_coefficient_long, aes(x = Comparison, y = simple_matching_coefficient_Index, color = Comparison)) +
    geom_boxplot(position = "dodge", alpha = 0.3, fill = "cornsilk",
                 outlier.shape = NA, size = 0.3) +
    geom_point(aes(shape = as.factor(Iteration)),
               position = position_jitter(width = 0.1), alpha = 0.7, size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3", linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", linewidth = 0.3, alpha = 0.5) +
    scale_color_colorblind() +
    scale_shape_manual(values = 0:19) +
    theme_light() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black", size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11),
      strip.text.y = element_text(angle = 90),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      x = "Comparison",
      y = "Simple matching coefficient",
      title = plot_title,
      subtitle = plot_subtitle
    ) +
    guides(shape = "none", color = "none")
  
  # Add annotations if requested
  if (add_annotations) {
    # Calculate statistics for each comparison type, correction type, AND parameter combination
    stats_by_group <- simple_matching_coefficient_long %>%
      group_by(Comparison, Correction_Type, ParameterCombination) %>%
      summarise(
        median_val = median(simple_matching_coefficient_Index, na.rm = TRUE),
        variance_val = var(simple_matching_coefficient_Index, na.rm = TRUE),
        min_val = min(simple_matching_coefficient_Index, na.rm = TRUE),
        mode_val = calculate_mode(simple_matching_coefficient_Index),
        .groups = "drop"
      )
    
    # Get unique comparisons in the order they appear in the plot
    unique_comparisons <- levels(simple_matching_coefficient_long$Comparison)
    if (is.null(unique_comparisons)) {
      unique_comparisons <- unique(simple_matching_coefficient_long$Comparison)
    }
    
    # Create annotation data frame for each comparison, correction type, and parameter combination
    annotation_df <- data.frame()
    
    for (i in seq_along(unique_comparisons)) {
      comp <- unique_comparisons[i]
      
      # Get stats for this comparison across correction types and parameter combinations
      comp_stats <- stats_by_group[stats_by_group$Comparison == comp,]
      
      if (nrow(comp_stats) > 0) {
        # Create annotations for this comparison for EACH correction type and parameter combination
        for (j in 1:nrow(comp_stats)) {
          correction_type <- comp_stats$Correction_Type[j]
          param_combo <- comp_stats$ParameterCombination[j]
          stats <- comp_stats[j,]
          
          comp_annotations <- data.frame(
            Comparison = rep(comp, 4),
            Correction_Type = rep(correction_type, 4),
            ParameterCombination = rep(param_combo, 4), # Include parameter combination
            x = c(i - 0.1, i + 0.1, i - 0.1, i + 0.1), # Side by side positioning
            y = c(1.5, 1.5, 1.15, 1.15), # Two rows
            color_group = rep(comp, 4),
            label = c(
              paste("Med:", round(stats$median_val, 3)),
              paste("Var:", round(stats$variance_val, 3)),
              paste("Min:", round(stats$min_val, 3)),
              paste("Mode:", round(stats$mode_val, 3))
            )
          )
          annotation_df <- rbind(annotation_df, comp_annotations)
        }
      }
    }
    
    # Add annotations to plot
    if (nrow(annotation_df) > 0) {
      p <- p + geom_text(
        data = annotation_df,
        aes(x = x, y = y, label = label, color = color_group),
        inherit.aes = FALSE,
        hjust = 1, vjust = 0.5, size = 2.2, angle = 90, fontface = "plain"
      )
    }
    
    # Adjust y-axis limits to accommodate annotations
    p <- p + scale_y_continuous(limits = c(0, 1.55), breaks = seq(0, 1, 0.25))
  } else {
    # Original y-axis limits without annotations
    p <- p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))
  }
  
  # Add faceting if requested and data supports it
  if (facet_by_param) {
    if (length(unique(simple_matching_coefficient_long$ParameterCombination)) > 1) {
      p <- p + facet_grid(Correction_Type ~ ParameterCombination, scales = "free_x") +
        labs(caption = "Rows show p-value correction type; Columns show parameter combinations")
    } else {
      # Single parameter combination - just facet by correction type
      p <- p + facet_wrap(Correction_Type ~ ., ncol = 1, strip.position = "right")
    }
  }
  
  return(p)
}

# ===============================================================================
# CUSTOM SINGLE EXPERIMENT FUNCTION
# ===============================================================================

# Additional function for custom single runs
run_custom_experiment <- function(CheckRemoved = FALSE, CheckThreefold = FALSE, OptimizeBetween = FALSE, NonNoiseSelection = FALSE,
                                  custom_nTrials = NULL) {

  if (is.null(custom_nTrials)) {
    custom_nTrials <- nTrials # Use global nTrials if not specified
  }

  cat(paste0("Running custom experiment with nTrials = ", custom_nTrials,
             " and parameters: CR=", CheckRemoved, ", CT=", CheckThreefold, ", OB=", OptimizeBetween, ", NNS=", NonNoiseSelection, "\n"))

  # Temporarily store original nTrials and set custom value
  original_nTrials <- nTrials
  assign("nTrials", custom_nTrials, envir = .GlobalEnv)

  # Run the experiment
  result <- run_experiment(CheckRemoved, CheckThreefold, OptimizeBetween, NonNoiseSelection)

  # Restore original nTrials
  assign("nTrials", original_nTrials, envir = .GlobalEnv)

  return(result)
}

# Usage:
# run_custom_experiment(custom_nTrials = 10, NonNoiseSelection = TRUE)

