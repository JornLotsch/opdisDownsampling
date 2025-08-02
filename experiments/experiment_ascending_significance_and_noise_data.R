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
nTrials <- 4 # Number of optimization trials per downsampling iteration
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
use_ABC_for_feature_selection <- TRUE # Whether to regard only the A subset features for the Jaccard Index

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
ascending_significance_and_noise_data <- generate_comprehensive_synthetic_data()

# Generic name for data
actual_data <- ascending_significance_and_noise_data
actual_class <- "Cls"

# Plot p-values from t-tests for each variable (may be slow for large datasets)
cat("Plotting p-values from t-tests...\n")
p_values <- apply(actual_data[, -1], 2, function(x) t.test(x ~ actual_data$Cls)$p.value)
plot(-log10(p_values), main = "-log10(p-values) from t-tests by variable", ylab = "p-value", xlab = "Variable index")

# Determine trial configurations to run
if (nTrials != 1) {
  trial_configs <- c(1, nTrials) # Run both single trial and full trials
} else {
  trial_configs <- nTrials # Only run single trial
}

# Execute experiments based on trial configuration
for (current_trial_count in trial_configs) {
  cat(sprintf("Running experiment with %d trial(s)...\n", current_trial_count))

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

# Print summary of results
if (exists("results_experiment_ascending_significance_1Trial")) {
  cat("Single trial results available in: results_experiment_ascending_significance_1Trial\n")
}
if (exists("results_experiment_ascending_significance_nTrials")) {
  cat(sprintf("Multiple trial results available in: results_experiment_ascending_significance_nTrials (%d combinations)\n",
              length(results_experiment_ascending_significance_nTrials)))
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
    # Calculate statistics for each comparison type
    stats_by_comparison <- correlation_data %>%
      group_by(Comparison) %>%
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

    # Create annotation data frame - 2x2 matrix arrangement for each comparison
    annotation_df <- data.frame()

    for (i in seq_along(unique_comparisons)) {
      comp <- unique_comparisons[i]
      stats <- stats_by_comparison[stats_by_comparison$Comparison == comp,]

      if (nrow(stats) > 0) {
        # Create annotations for this comparison (side by side positioning)
        comp_annotations <- data.frame(
          Comparison = rep(comp, 4),
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

    # Add annotations to plot
    if (nrow(annotation_df) > 0) {
      p <- p + geom_text(
        data = annotation_df,
        aes(x = x, y = y, label = label, color = color_group),
        inherit.aes = FALSE,
        hjust = 1, vjust = 0.5, size = 2.5, angle = 90, fontface = "plain"
      )
    }
  } else {
    # Original y-axis limits without annotations
    p <- p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))
  }

  # Add faceting if requested and multiple parameter combinations exist
  if (facet_by_param && length(unique(correlation_data$ParameterCombination)) > 1) {
    p <- p + facet_wrap(. ~ ParameterCombination, nrow = 1) +
      guides(shape = guide_legend(nrow = 1), color = "none")
  } else {
    p <- p + guides(color = "none", shape = "none")
  }

  if (add_annotations) {
    scale_y_continuous(
      limits = c(floor(min(correlation_data$Tau, na.rm = TRUE) * 4) / 4, 1.55), # Extended upper limit for annotations
      breaks = seq(from = floor(min(correlation_data$Tau, na.rm = TRUE) * 4) / 4,
                   to = 1, by = 0.25)
      )
  } else {
    scale_y_continuous(
        limits = c(min(correlation_data$Tau), 1), # Extended upper limit for annotations
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
create_jaccard_plot <- function(jaccard_data, plot_title = "Jaccard Index Analysis",
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
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3", size = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", size = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 0.3, alpha = 0.5) +
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
      y = "Jaccard Index",
      title = plot_title,
      subtitle = plot_subtitle
    ) +
    guides(shape = "none", color = "none")

  # Add annotations if requested
  if (add_annotations) {
    # Calculate statistics for each comparison type and correction type
    stats_by_group <- jaccard_long %>%
      group_by(Comparison, Correction_Type) %>%
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

    # Create annotation data frame - 2x2 matrix arrangement for each comparison
    annotation_df <- data.frame()

    for (i in seq_along(unique_comparisons)) {
      comp <- unique_comparisons[i]

      # Get stats for this comparison across correction types
      comp_stats <- stats_by_group[stats_by_group$Comparison == comp,]

      if (nrow(comp_stats) > 0) {
        # Create annotations for this comparison (side by side positioning)
        comp_annotations <- data.frame(
          Comparison = rep(comp, 4 * nrow(comp_stats)),
          Correction_Type = rep(comp_stats$Correction_Type, each = 4),
          x = rep(c(i - 0.1, i + 0.1, i - 0.1, i + 0.1), nrow(comp_stats)), # Side by side positioning
          y = rep(c(1.5, 1.5, 1.15, 1.15), nrow(comp_stats)), # Two rows
          color_group = rep(comp, 4 * nrow(comp_stats)),
          label = c(
            sapply(1:nrow(comp_stats), function(j) {
              c(
                paste("Med:", round(comp_stats$median_val[j], 3)),
                paste("Var:", round(comp_stats$variance_val[j], 3)),
                paste("Min:", round(comp_stats$min_val[j], 3)),
                paste("Mode:", round(comp_stats$mode_val[j], 3))
              )
            })
          )
        )
        annotation_df <- rbind(annotation_df, comp_annotations)
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
      p <- p + facet_wrap(~Correction_Type, ncol = 1)
    }
  }

  return(p)
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
  single_correlations <- all_correlations[all_correlations$ParameterCombination == "OpF_NoF",]
  multiple_correlations <- all_correlations

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
  single_jaccard <- all_jaccard_data[all_jaccard_data$ParameterCombination == "OpF_NoF",]
  multiple_jaccard <- all_jaccard_data

  p_single_trial_jaccard <- create_jaccard_plot(
    single_jaccard,
    plot_title = "Single Trial",
    plot_subtitle = paste0("1 trial per iteration,\n", nSamples, " iterations"),
    facet_by_param = FALSE,
    add_annotations = add_annotations
  )

  p_multiple_jaccard <- create_jaccard_plot(
    multiple_jaccard,
    plot_title = "Feature Selection Overlap by P-value Correction Method",
    plot_subtitle = paste("Analysis with", format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = TRUE,
    add_annotations = add_annotations
  )

} else if (nrow(all_jaccard_data) > 0) {
  # Create single comprehensive plot
  p_all_jaccard <- create_jaccard_plot(
    all_jaccard_data,
    plot_title = "Jaccard Index Analysis Results",
    plot_subtitle = paste("Analysis with", format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    facet_by_param = has_multiple_params,
    add_annotations = add_annotations
  )
}

# Combine plots for comprehensive visualization
if (exists("p_single_trial_correlations") && exists("p_multiple_correlations") &&
    exists("p_single_trial_jaccard") && exists("p_multiple_jaccard")) {

  p_combined_analysis <- cowplot::plot_grid(
    p_single_trial_correlations,
    p_multiple_correlations,
    p_single_trial_jaccard,
    p_multiple_jaccard,
    labels = "AUTO",
    nrow = 2,
    rel_widths = c(1, 4),
    align = "hv", axis = "tb"
  )

} else if (exists("p_all_correlations") && exists("p_all_jaccard")) {

  p_combined_analysis <- cowplot::plot_grid(
    p_all_correlations,
    p_all_jaccard,
    labels = c("A", "B"),
    nrow = 2,
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
  base_width <- 10 # Width for single parameter combination
  width_per_column <- 5 # Additional width per extra parameter combination
  plot_width <- base_width + (n_param_combinations - 1) * width_per_column

  # Set reasonable bounds
  plot_width <- max(8, min(plot_width, 24)) # Between 8 and 24 inches

  # Save with descriptive filename
  output_filename <- paste0("combined_analysis_comparison_",
                            format(nTrials, scientific = FALSE), "trials_",
                            nSamples, "iterations_", downsampling_size, "sampled.svg")

  ggsave(output_filename, p_combined_analysis, width = plot_width, height = 18)
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

# Save all 
save.image(file = "experiment_ascending_significance_and_noise_data.RData")
