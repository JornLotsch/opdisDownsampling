# ===============================================================================
# ARTIFICIAL DATA SETS GENERATION AND EVALUATION SCRIPT
# ===============================================================================
#
# Purpose: Performs assessment of data splitting parameter effects on feature
#          significance preservation.
#
# ===============================================================================

# ===============================================================================
# EXPERIMENTAL CONFIGURATION PARAMETERS
# ===============================================================================
random_seed <- 42 # Random seed for reproducibility
significance_threshold <- 0.05 # Statistical significance threshold for p-values
downsampling_size <- 0.8 # Proportion of original data to retain (80%)
nSamples <- 10 # Number of downsampling iterations to perform
nTrials <- 10000 # Number of optimization trials per downsampling iteration
TestStat <- "ad" # Statistical test for distribution comparison (Anderson-Darling)

# opdisDownsampling-specific parameters
UniformTestStat <- "ks" # Statistical test for non-uniform variable selection (Kolmogorov-Smirnov)
UniformThreshold <- 0.05 # Threshold value for non-uniform variable selection
JobSize <- 0 # Number of seeds per chunk (0 = no chunking)
verbose <- FALSE # Whether to print chunk size diagnostics
WorstSample <- FALSE # Test if results worsen when reversing selection criterion
PCAimportance <- FALSE # Whether to use PCA to identify relevant variables

# Analysis and visualization parameters
use_ABC_for_Pvalues <- TRUE # Calculate Tau only for relevant features
use_ABC_for_feature_selection <- TRUE # Whether to regard only the A subset features for the Jaccard index

# Default values for data splitting optimization options
CheckRemoved <- FALSE # Whether to optimize for reduced AND removed sets
CheckThreefold <- FALSE # Whether to optimize for reduced AND removed AND between both sets
OptimizeBetween <- FALSE # Whether to optimize between reduced and removed sets
NonNoiseSelection <- FALSE # Whether to use noise detection for variable pre-selection
use_regression_for_p <- FALSE # Whether to do statistics for variable importance using U tests or regression
add_annotations <- TRUE # Whether to annotate plots with summary stats of results
split_method <- "opdis" # Which split method to use (currently only "opdis")

# Execution control
RUN_PARAMETER_IMPACT_ANALYSIS <- FALSE # Set to FALSE to skip parameter impact analysis

# Directory paths
EXPERIMENTS_DIR <- "/home/joern/Aktuell/DownSamplingStructure/12RLibrary/opdisDownsampling/experiments"

# ===============================================================================
# SET PATH AND LOAD FUNCTIONS
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

# Check if required source files exist before sourcing
required_files <- c("generate_ascending_significance_and_noise_data.R",
                    "experiments_functions.R")

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required files not found in current directory (", getwd(), "):\n  - ",
       paste(missing_files, collapse = "\n  - "),
       "\n\nPlease ensure these files are in the working directory or set EXPERIMENTS_DIR environment variable.")
}

cat("All required files found. Sourcing...\n")

# Source functions with better error handling
for (file in required_files) {
  tryCatch({
    source(file)
    cat("Successfully sourced:", file, "\n")
  }, error = function(e) {
    stop("Error sourcing '", file, "': ", e$message)
  })
}

cat("All files sourced successfully.\n")


# ===============================================================================
# EXPERIMENT EXECUTION LOGIC
# ===============================================================================

# Generate all parameter combinations for comprehensive testing
params <- expand.grid(
# CheckRemoved = c(TRUE, FALSE), # Test both optimization strategies
# CheckRemoved  = c(TRUE, FALSE), # Test both optimization strategies
  OptimizeBetween = c(TRUE, FALSE), # Test both optimization strategies
  NonNoiseSelection = c(TRUE, FALSE) # Test both variable selection approaches
)

# Input validation for number of trials
if (!is.numeric(nTrials) || nTrials <= 0) {
  stop("nTrials must be a positive integer")
}

# Define all possible opdisDownsampling parameters and their defaults
all_possible_params <- list(
  CheckRemoved = FALSE,
  CheckThreefold = FALSE,
  OptimizeBetween = FALSE,
  NonNoiseSelection = FALSE,
  PCAimportance = FALSE,
  UniformTestStat = "ks",
  UniformThreshold = 0.05,
  JobSize = 0,
  verbose = FALSE,
  WorstSample = FALSE
)

# Function to create parameter list with defaults for missing parameters
create_full_param_list <- function(param_row, all_defaults) {
  # Start with all defaults
  full_params <- all_defaults

  # Override with values from params data frame
  param_names <- names(param_row)
  for (name in param_names) {
    if (name %in% names(full_params)) {
      full_params[[name]] <- param_row[[name]]
    }
  }

  return(full_params)
}

# Generate the synthetic data
# ascending_significance_and_noise_data <- generate_custom_biomedical_data(n_normal = 40, n_uniform = 40)
ascending_significance_and_noise_data <- generate_biomedical_synthetic_data(num_variables = 15)

# Generic name for data
actual_data <- ascending_significance_and_noise_data
actual_class <- "Cls"

# Plot p-values from t-tests for each variable (may be slow for large datasets)
cat("Plotting p-values from t-tests...\n")
p_values <- apply(actual_data[, -1], 2, function(x) t.test(x ~ actual_data$Cls)$p.value)
plot(-log10(p_values), main = "-log10(p-values) from t-tests by variable", ylab = "p-value", xlab = "Variable index")
# par(mfrow=c(2,1))
# barplot(-log10(p_values), main = "-log10(p-values) from t-tests by variable", ylab = "p-value", xlab = "Variable index", las = 2)
# par(mfrow=c(1,1))


# Determine trial configurations to run
if (nTrials != 1) {
  trial_configs <- c(1, nTrials) # Run both single trial and full trials
} else {
  trial_configs <- nTrials # Only run single trial
}

# Execute experiments based on trial configuration
for (current_trial_count in trial_configs) {
  cat(sprintf("Running experiment with %d trial(s)...\n", current_trial_count))
  nTrials = current_trial_count

  if (current_trial_count == 1) {
    # Single run with all parameters set to FALSE - quick test mode
    cat("Running single trial with all parameters set to FALSE...\n")

    # Create full parameter list with all defaults (all FALSE)
    single_trial_params <- all_possible_params

    results_experiment_ascending_significance_1Trial <- run_experiment(
      data_df = actual_data,
      class_name = actual_class,
      CheckRemoved = single_trial_params$CheckRemoved,
      CheckThreefold = single_trial_params$CheckThreefold,
      OptimizeBetween = single_trial_params$OptimizeBetween,
      NonNoiseSelection = single_trial_params$NonNoiseSelection,
      PCAimportance = single_trial_params$PCAimportance,
      UniformTestStat = single_trial_params$UniformTestStat,
      UniformThreshold = single_trial_params$UniformThreshold,
      JobSize = single_trial_params$JobSize,
      verbose = single_trial_params$verbose,
      WorstSample = single_trial_params$WorstSample
    )

  } else if (current_trial_count > 1) {
    # Run all parameter combinations - comprehensive analysis mode
    cat(sprintf("Running comprehensive analysis with %d parameter combinations...\n", nrow(params)))

    # Print parameter combinations that will be tested
    cat("Parameter combinations to test:\n")
    for (i in seq_len(nrow(params))) {
      param_summary <- paste(names(params), "=", params[i,], collapse = ", ")
      cat(sprintf("  %d: %s\n", i, param_summary))
    }

    results_experiment_ascending_significance_nTrials <-
      lapply(seq_len(nrow(params)), function(i) {
        # Create full parameter list for this combination
        current_params <- create_full_param_list(params[i,], all_possible_params)

        # Create parameter summary for logging
        varied_params <- names(params)
        param_values <- paste(varied_params, "=", params[i, varied_params], collapse = ", ")

        cat(sprintf("  Running combination %d/%d: %s\n", i, nrow(params), param_values))

        res <- run_experiment(
          data_df = actual_data,
          class_name = actual_class,
          CheckRemoved = current_params$CheckRemoved,
          CheckThreefold = current_params$CheckThreefold,
          OptimizeBetween = current_params$OptimizeBetween,
          NonNoiseSelection = current_params$NonNoiseSelection,
          PCAimportance = current_params$PCAimportance,
          UniformTestStat = current_params$UniformTestStat,
          UniformThreshold = current_params$UniformThreshold,
          JobSize = current_params$JobSize,
          verbose = current_params$verbose,
          WorstSample = current_params$WorstSample
        )
        return(res)
      })

    # Create names for results based on the parameters that vary in the params data frame
    param_names <- names(params)
    result_names <- apply(params, 1, function(row) {
      param_strings <- paste0(
        substr(param_names, 1, 2), # Use first 2 letters of parameter name
        ifelse(row, "T", "F") # T for TRUE, F for FALSE
      )
      paste(param_strings, collapse = "_")
    })

    names(results_experiment_ascending_significance_nTrials) <- result_names

    cat("Parameter combinations completed:\n")
    for (i in seq_along(names(results_experiment_ascending_significance_nTrials))) {
      cat(sprintf("  %s\n", names(results_experiment_ascending_significance_nTrials)[i]))
    }
  }
}

cat("All analyses complete!\n")

# Print summary of results
if (exists("results_experiment_ascending_significance_1Trial")) {
  cat("Single trial results available in: results_experiment_ascending_significance_1Trial\n")
}
if (exists("results_experiment_ascending_significance_nTrials")) {
  cat(sprintf("Multiple trial results available in: results_experiment_ascending_significance_nTrials (%d combinations)\n",
              length(results_experiment_ascending_significance_nTrials)))
}


# ---------------------------------------------------------------------------
# MAIN POST-PROCESSING EXECUTION
# ---------------------------------------------------------------------------

cat("Starting post-processing and comparative analysis...\n")

# Extract correlation data from available results
single_trial_results <- if (exists("results_experiment_ascending_significance_1Trial")) {
  results_experiment_ascending_significance_1Trial
} else {
  NULL
}

multiple_trial_results <- if (exists("results_experiment_ascending_significance_nTrials")) {
  results_experiment_ascending_significance_nTrials
} else {
  NULL
}

# Extract and combine all correlation data
all_correlations <- extract_correlation_data(
  results_single = single_trial_results,
  results_multiple = multiple_trial_results
)

# Extract and combine all Jaccard data
all_jaccard_data <- extract_jaccard_data(
  results_single = single_trial_results,
  results_multiple = multiple_trial_results
)

# Extract and combine all simple matching coefficient data
all_simple_matching_coefficient_data <- extract_simple_matching_coefficient_data(
  results_single = single_trial_results,
  results_multiple = multiple_trial_results
)

# ===============================================================================
# COMPARATIVE VISUALIZATION
# ===============================================================================

cat("Creating comparative visualizations...\n")

# Determine if we have multiple parameter combinations for appropriate plotting
has_multiple_params <- length(unique(all_correlations$ParameterCombination)) > 1
has_single_trial <- !is.null(single_trial_results)
has_multiple_trials <- !is.null(multiple_trial_results)

# Create correlation plots
if (has_single_trial && has_multiple_trials) {
  # Create separate plots for single and multiple trials
  single_correlations <- all_correlations[all_correlations$SingleMultiple == "Single",]
  multiple_correlations <- all_correlations[all_correlations$SingleMultiple == "Multiple",]

  p_single_trial_correlations <- create_correlation_plot(
    single_correlations,
    plot_title = "Single Trial",
    plot_subtitle = paste0("1 trial per iteration,\n", nSamples, " iterations"),
    facet_by_param = FALSE,
    add_annotations = add_annotations
  )

  p_multiple_correlations <- create_correlation_plot(
    multiple_correlations,
    plot_title = "Multiple Parameter Combinations",
    plot_subtitle = paste(format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = TRUE,
    add_annotations = add_annotations
  )

} else if (nrow(all_correlations) > 0) {
  # Create single comprehensive plot
  p_all_correlations <- create_correlation_plot(
    all_correlations,
    plot_title = "Correlation Analysis Results",
    plot_subtitle = paste("Analysis with", format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = has_multiple_params,
    add_annotations = add_annotations
  )
}

# Create Jaccard plots
if (has_single_trial && has_multiple_trials) {
  # Create separate plots for single and multiple trials
  single_jaccard <- all_jaccard_data[all_jaccard_data$SingleMultiple == "Single",]
  multiple_jaccard <- all_jaccard_data[all_jaccard_data$SingleMultiple == "Multiple",]

  p_single_trial_jaccard <- create_jaccard_plot(
    single_jaccard,
    plot_title = "Single trial",
    plot_subtitle = paste0("1 trial per iteration,\n", nSamples, " iterations"),
    facet_by_param = TRUE,
    add_annotations = add_annotations
  )

  p_multiple_jaccard <- create_jaccard_plot(
    multiple_jaccard,
    plot_title = "Feature selection overlap by p-value correction method",
    plot_subtitle = paste("Analysis with", format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = TRUE,
    add_annotations = add_annotations
  )

} else if (nrow(all_jaccard_data) > 0) {
  # Create single combined plot
  p_all_jaccard <- create_jaccard_plot(
    all_jaccard_data,
    plot_title = "Jaccard index analysis results",
    plot_subtitle = paste("Analysis with", format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = has_multiple_params,
    add_annotations = add_annotations
  )
}

# Create simple matching coefficient plots
if (has_single_trial && has_multiple_trials) {
  # Create separate plots for single and multiple trials
  single_simple_matching_coefficient <- all_simple_matching_coefficient_data[all_simple_matching_coefficient_data$SingleMultiple == "Single",]
  multiple_simple_matching_coefficient <- all_simple_matching_coefficient_data[all_simple_matching_coefficient_data$SingleMultiple == "Multiple",]

  p_single_trial_simple_matching_coefficient <- create_simple_matching_coefficient_plot(
    single_simple_matching_coefficient,
    plot_title = "Single trial",
    plot_subtitle = paste0("1 trial per iteration,\n", nSamples, " iterations"),
    facet_by_param = TRUE,
    add_annotations = add_annotations
  )

  p_multiple_simple_matching_coefficient <- create_simple_matching_coefficient_plot(
    multiple_simple_matching_coefficient,
    plot_title = "Feature selection overlap by p-value correction method",
    plot_subtitle = paste("Analysis with", format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = TRUE,
    add_annotations = add_annotations
  )

} else if (nrow(all_simple_matching_coefficient_data) > 0) {
  # Create single combined plot
  p_all_simple_matching_coefficient <- create_simple_matching_coefficient_plot(
    all_simple_matching_coefficient_data,
    plot_title = "Simple matching coefficient analysis results",
    plot_subtitle = paste("Analysis with", format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = has_multiple_params,
    add_annotations = add_annotations
  )
}

# Combine plots for overview
if (exists("p_single_trial_correlations") && exists("p_multiple_correlations") &&
    exists("p_single_trial_jaccard") && exists("p_multiple_jaccard") &&
    exists("p_single_trial_simple_matching_coefficient") && exists("p_multiple_simple_matching_coefficient")) {

  p_combined_analysis <- cowplot::plot_grid(
    p_single_trial_correlations,
    p_multiple_correlations,
    p_single_trial_jaccard,
    p_multiple_jaccard,
    p_single_trial_simple_matching_coefficient,
    p_multiple_simple_matching_coefficient,
    labels = "AUTO",
    nrow = 3,
    rel_widths = c(1, 4),
    align = "hv", axis = "tb"
  )

} else if (exists("p_all_correlations") && exists("p_all_jaccard") && exists("p_all_simple_matching_coefficient")) {

  p_combined_analysis <- cowplot::plot_grid(
    p_all_correlations,
    p_all_jaccard,
    p_all_simple_matching_coefficient,
    labels = "AUTO",
    nrow = 3,
    align = "hv", axis = "tb"
  )
}

# Display and save results
if (exists("p_combined_analysis")) {
  print(p_combined_analysis)

  # Calculate dynamic width based on number of facets in p_multiple_correlations
  # Extract number of parameter combinations from the correlation data
  if (exists("combined_correlations") && nrow(combined_correlations) > 0) {
    n_param_combinations <- length(unique(combined_correlations$ParameterCombination))
  } else {
    n_param_combinations <- 1 # fallback
  }

  # Calculate width: base width + extra width per additional column
  base_width <- 12 # Width for single parameter combination
  width_per_column <- 5 # Additional width per extra parameter combination
  plot_width <- base_width + (n_param_combinations - 1) * width_per_column

  # Set reasonable bounds
  plot_width <- max(8, min(plot_width, 24)) # Between 8 and 24 inches

  # Save with descriptive filename
  output_filename <- paste0("combined_analysis_comparison_",
                            format(nTrials, scientific = FALSE), "trials_",
                            nSamples, "iterations_", downsampling_size, "sampled.svg")

  ggsave(output_filename, p_combined_analysis, width = plot_width, height = 25)
  cat(sprintf("Combined analysis plot saved as: %s (width: %.1f inches, %d param combinations)\n",
              output_filename, plot_width, n_param_combinations))
}

# Print summary of processed results
cat("\n=== POST-PROCESSING SUMMARY ===\n")
if (nrow(all_correlations) > 0) {
  cat(sprintf("Processed correlation data: %d observations from %d parameter combination(s)\n",
              nrow(all_correlations), length(unique(all_correlations$ParameterCombination))))
}
if (nrow(all_jaccard_data) > 0) {
  cat(sprintf("Processed Jaccard data: %d observations from %d parameter combination(s)\n",
              nrow(all_jaccard_data), length(unique(all_jaccard_data$ParameterCombination))))
}

if (exists("results_experiment_ascending_significance_1Trial")) {
  cat("Single trial results processed: results_experiment_ascending_significance_1Trial\n")
}
if (exists("results_experiment_ascending_significance_nTrials")) {
  cat(sprintf("Multiple trial results processed: results_experiment_ascending_significance_nTrials (%d combinations)\n",
              length(results_experiment_ascending_significance_nTrials)))
}

cat("Post-processing and comparative analysis complete!\n")



if (RUN_PARAMETER_IMPACT_ANALYSIS) {
  # ===============================================================================
  # PARAMETER IMPACT ANALYSIS
  # ===============================================================================

  # Check for identical parameter combinations in all_correlations
  check_identical_parameter_combinations <- function(df) {
    # Get unique parameter combinations
    param_combinations <- unique(df$ParameterCombination)

    # Extract Tau values for each parameter combination
    tau_blocks <- lapply(param_combinations, function(param) {
      subset_df <- df[df$ParameterCombination == param,]
      # Sort by iteration to ensure correct order
      subset_df <- subset_df[order(subset_df$Iteration),]
      return(subset_df$Tau)
    })
    names(tau_blocks) <- param_combinations

    # Find identical blocks
    identical_groups <- list()
    processed <- character(0)

    for (i in 1:length(tau_blocks)) {
      param1 <- names(tau_blocks)[i]

      # Skip if already processed
      if (param1 %in% processed) next

      # Find all parameters with identical Tau sequences
      identical_params <- c(param1)

      for (j in (i + 1):length(tau_blocks)) {
        if (j > length(tau_blocks)) break

        param2 <- names(tau_blocks)[j]

        # Skip if already processed
        if (param2 %in% processed) next

        # Check if Tau sequences are identical
        if (identical(tau_blocks[[param1]], tau_blocks[[param2]])) {
          identical_params <- c(identical_params, param2)
          processed <- c(processed, param2)
        }
      }

      # Add to groups if more than one parameter has identical sequences
      if (length(identical_params) > 1) {
        identical_groups[[length(identical_groups) + 1]] <- identical_params
      }

      processed <- c(processed, param1)
    }

    return(list(
      tau_blocks = tau_blocks,
      identical_groups = identical_groups
    ))
  }

  # Check for identical parameter combinations
  identical_check <- check_identical_parameter_combinations(all_correlations)

  # Print results
  cat("\n=== PARAMETER IMPACT ANALYSIS ===\n")
  cat("Identical Parameter Combination Groups:\n")
  if (length(identical_check$identical_groups) > 0) {
    for (i in 1:length(identical_check$identical_groups)) {
      cat(paste("Group", i, ":\n"))
      for (param in identical_check$identical_groups[[i]]) {
        cat(paste("  -", param, "\n"))
      }
      cat("\n")
    }
  } else {
    cat("No identical parameter combinations found.\n")
  }

  # Show first few Tau values for each parameter combination for verification
  cat("\nFirst 3 Tau values for each parameter combination:\n")
  for (param in names(identical_check$tau_blocks)) {
    tau_values <- identical_check$tau_blocks[[param]]
    cat(paste(param, ":", paste(round(tau_values[1:min(3, length(tau_values))], 4), collapse = ", "), "...\n"))
  }

  # Analyze identical groups for best performance
  analyze_identical_groups <- function(identical_check) {
    if (length(identical_check$identical_groups) == 0) {
      cat("No identical groups found.\n")
      return(NULL)
    }

    group_stats <- data.frame(
      Group = integer(),
      Parameters = character(),
      MedianTau = numeric(),
      MinTau = numeric(),
      MaxTau = numeric(),
      MeanTau = numeric(),
      stringsAsFactors = FALSE
    )

    cat("Analysis of identical parameter combination groups:\n\n")

    for (i in 1:length(identical_check$identical_groups)) {
      group_params <- identical_check$identical_groups[[i]]

      # Get Tau values for this group (they're all identical, so use first parameter)
      tau_values <- identical_check$tau_blocks[[group_params[1]]]

      # Calculate statistics
      median_tau <- median(tau_values)
      min_tau <- min(tau_values)
      max_tau <- max(tau_values)
      mean_tau <- mean(tau_values)

      # Store in data frame
      group_stats <- rbind(group_stats, data.frame(
        Group = i,
        Parameters = paste(group_params, collapse = "; "),
        MedianTau = median_tau,
        MinTau = min_tau,
        MaxTau = max_tau,
        MeanTau = mean_tau,
        stringsAsFactors = FALSE
      ))

      # Print group details
      cat(paste("Group", i, ":\n"))
      cat(paste("  Parameters:", paste(group_params, collapse = ", "), "\n"))
      cat(paste("  Median Tau:", round(median_tau, 4), "\n"))
      cat(paste("  Min Tau:", round(min_tau, 4), "\n"))
      cat(paste("  Max Tau:", round(max_tau, 4), "\n"))
      cat(paste("  Mean Tau:", round(mean_tau, 4), "\n"))
      cat(paste("  Tau values:", paste(round(tau_values, 4), collapse = ", "), "\n"))
      cat("\n")
    }

    # Find best groups
    best_median_group <- which.max(group_stats$MedianTau)
    best_min_group <- which.max(group_stats$MinTau)

    cat("SUMMARY:\n")
    cat("========\n")
    cat(paste("Group with HIGHEST MEDIAN Tau:", best_median_group,
              "(Median =", round(group_stats$MedianTau[best_median_group], 4), ")\n"))
    cat(paste("  Parameters:", group_stats$Parameters[best_median_group], "\n"))
    cat("\n")

    cat(paste("Group with HIGHEST MINIMUM Tau:", best_min_group,
              "(Min =", round(group_stats$MinTau[best_min_group], 4), ")\n"))
    cat(paste("  Parameters:", group_stats$Parameters[best_min_group], "\n"))

    if (best_median_group == best_min_group) {
      cat("\n*** The same group has both highest median and highest minimum Tau! ***\n")
    }

    return(group_stats)
  }

  # Run the analysis
  group_analysis <- analyze_identical_groups(identical_check)

  # Show the summary table
  if (!is.null(group_analysis)) {
    cat("\n\nSummary Table:\n")
    print(group_analysis)
  }

  # Identify redundant parameters based on identical results
  identify_redundant_parameters <- function(identical_check) {
    if (length(identical_check$identical_groups) == 0) {
      cat("No identical groups found - no redundant parameters to identify.\n")
      return(NULL)
    }

    # Parse parameter combinations to understand structure
    parse_parameter_combination <- function(param_combo) {
      # Extract individual parameters from string like "OpF_NoF"
      parts <- strsplit(param_combo, "_")[[1]]
      params <- list()

      for (part in parts) {
        if (grepl("^OB", part)) {
          params$OptimizeBetween <- grepl("TRUE", part)
        } else if (grepl("^NNS", part)) {
          params$NonNoiseSelection <- grepl("TRUE", part)
        }
      }

      return(params)
    }

    # Get all unique parameter combinations
    all_params <- unique(names(identical_check$tau_blocks))

    # Parse all parameter combinations
    parsed_params <- lapply(all_params, parse_parameter_combination)
    names(parsed_params) <- all_params

    # Analyze redundancy patterns
    redundancy_analysis <- list()

    cat("REDUNDANCY ANALYSIS:\n")
    cat("===================\n\n")

    for (i in 1:length(identical_check$identical_groups)) {
      group_params <- identical_check$identical_groups[[i]]

      cat(paste("Group", i, "- Identical Results:\n"))
      for (param in group_params) {
        cat(paste("  ", param, "\n"))
      }

      # Parse parameters for this group
      group_parsed <- parsed_params[group_params]

      # Find what parameters vary within this identical group
      param_names <- c("OptimizeBetween", "NonNoiseSelection")
      varying_params <- character()
      constant_params <- character()

      for (param_name in param_names) {
        values <- sapply(group_parsed, function(x) x[[param_name]])
        if (length(unique(values)) > 1) {
          varying_params <- c(varying_params, param_name)
        } else {
          constant_params <- c(constant_params, paste0(param_name, "=", values[1]))
        }
      }

      cat(paste("  Constant parameters:", paste(constant_params, collapse = ", "), "\n"))
      cat(paste("  Varying parameters:", paste(varying_params, collapse = ", "), "\n"))

      # Determine redundancy implications
      if (length(varying_params) == 1) {
        cat(paste("  -> REDUNDANCY FOUND:", varying_params[1], "has no effect when", paste(constant_params, collapse = ", "), "\n"))
      } else if (length(varying_params) > 1) {
        cat(paste("  -> INTERACTION:", paste(varying_params, collapse = " + "), "combinations have no effect when", paste(constant_params, collapse = ", "), "\n"))
      }

      cat("\n")

      # Store for summary
      redundancy_analysis[[i]] <- list(
        group = i,
        parameters = group_params,
        varying = varying_params,
        constant = constant_params,
        redundancy_type = if (length(varying_params) == 1) "single_parameter" else "interaction"
      )
    }

    # Summary of redundant parameters
    cat("REDUNDANCY SUMMARY:\n")
    cat("==================\n")

    redundant_params <- character()

    for (analysis in redundancy_analysis) {
      if (analysis$redundancy_type == "single_parameter") {
        redundant_param <- analysis$varying[1]
        context <- paste(analysis$constant, collapse = ", ")

        cat(paste("• Parameter '", redundant_param, "' is REDUNDANT when ", context, "\n", sep = ""))
        cat(paste("  (Switching ", redundant_param, " TRUE/FALSE produces identical results)\n", sep = ""))

        redundant_params <- c(redundant_params, redundant_param)

      } else if (analysis$redundancy_type == "interaction") {
        context <- paste(analysis$constant, collapse = ", ")
        varying_combo <- paste(analysis$varying, collapse = " + ")

        cat(paste("• Parameter interaction '", varying_combo, "' is redundant when ", context, "\n", sep = ""))
        cat(paste("  (Different combinations of ", varying_combo, " produce identical results)\n", sep = ""))
      }
    }

    # Find globally redundant parameters (redundant in all contexts)
    if (length(redundant_params) > 0) {
      param_counts <- table(redundant_params)
      total_groups <- length(identical_check$identical_groups)

      globally_redundant <- names(param_counts)[param_counts == total_groups]

      if (length(globally_redundant) > 0) {
        cat("\nGLOBALLY REDUNDANT PARAMETERS:\n")
        cat("==============================\n")
        for (param in globally_redundant) {
          cat(paste("• '", param, "' can be REMOVED - it has no effect in any tested combination\n", sep = ""))
        }
      }
    }

    # Recommendation
    cat("\nRECOMMENDATIONS:\n")
    cat("===============\n")

    if (length(redundant_params) > 0) {
      unique_redundant <- unique(redundant_params)
      cat("Parameters that can potentially be simplified or removed:\n")
      for (param in unique_redundant) {
        # Count in how many contexts this parameter is redundant
        contexts <- sum(sapply(redundancy_analysis, function(x) param %in% x$varying && x$redundancy_type == "single_parameter"))
        cat(paste("  - ", param, " (redundant in ", contexts, " context(s))\n", sep = ""))
      }
    } else {
      cat("No single parameters found to be completely redundant.\n")
      cat("However, some parameter combinations produce identical results.\n")
    }

    return(list(
      analysis = redundancy_analysis,
      redundant_parameters = unique(redundant_params),
      recommendations = "Check the printed analysis above"
    ))
  }

  # Run the redundancy analysis
  redundancy_results <- identify_redundant_parameters(identical_check)

  # Store results for later access
  assign("parameter_impact_analysis", list(
    identical_check = identical_check,
    group_analysis = group_analysis,
    redundancy_results = redundancy_results
  ), envir = .GlobalEnv)

  cat("\nParameter impact analysis results stored in 'parameter_impact_analysis' variable.\n")

} else {
  cat("\n=== PARAMETER IMPACT ANALYSIS SKIPPED ===\n")
  cat("Set RUN_PARAMETER_IMPACT_ANALYSIS = TRUE to enable this analysis.\n")
}

# Final completion message
cat("\n=== ALL PROCESSING COMPLETE ===\n")
cat("Results stored in workspace variables for further analysis.\n")

# Only print parameter_impact_analysis if it was run
if (RUN_PARAMETER_IMPACT_ANALYSIS && exists("parameter_impact_analysis")) {
  print(parameter_impact_analysis)
}

# Analyze results per varibale class
source("analyse_variable_classes.R")

# Save all 
save.image(file = "experiment_ascending_significance_and_noise_data.RData")
