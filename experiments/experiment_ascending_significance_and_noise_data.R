#' Artificial Data Sets Generation and Evaluation Script
#'
#' Generates synthetic datasets for evaluating data splitting methods,
#' featuring variables with clear class separation and added noise variables.
#' Performs assessment of data splitting parameter effects on feature significance preservation 

# Attempt to set working directory to script location (for RStudio)
tryCatch({
  if (exists("rstudioapi::getSourceEditorContext", mode = "function")) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}, error = function(e) {
  message("Unable to set working directory automatically. Please set it manually if needed.")
})

setwd("/home/joern/Aktuell/DownSamplingStructure/12RLibrary/opdisDownsampling/experiments")

# Load required libraries
library(stats)      # Statistical functions
library(ggplot2)    # Data visualization
library(ggthemes)   # Additional ggplot2 themes
library(dplyr)      # Data manipulation
library(parallel)      # Data manipulation
library(pbmcapply)      # Data manipulation
library(cowplot)      # Data manipulation

# Configuration parameters
random_seed <- 42            # Random seed for reproducibility
samples_per_class <- 100     # Samples per class
num_variables <- 50          # Number of variables (features)
significance_threshold <- 0.05
downsampling_size <- 0.8
nSamples <- 10               # Number of downsampling iterations
nTrials <- 100000             # Number of trials per downsampling iteration
TestStat = "ad"

# Default values for data splitting options
CheckRemoved <- FALSE
OptimizeBetween <- FALSE
NonNoiseSelection <- FALSE

# Function to run experiment given parameters
run_experiment <- function(CheckRemoved, OptimizeBetween, NonNoiseSelection) {
  
  # Print current parameter combination
  print(c(CheckRemoved, OptimizeBetween, NonNoiseSelection))
  
  #' ===== Dataset Generation Functions =====
  
  #' Generate synthetic dataset with ascending class separation
  generate_synthetic_data <- function(n_samples, n_vars, seed) {
    class1_means <- seq(1, 50, length.out = n_vars)
    class2_means <- seq(1, 48, length.out = n_vars)
    
    variables <- lapply(1:n_vars, function(i) {
      set.seed(seed)
      class1_data <- rnorm(n_samples, mean = class1_means[i], sd = 5)
      class2_data <- rnorm(n_samples, mean = class2_means[i], sd = 5)
      c(class1_data, class2_data)
    })
    
    data.frame(
      Cls = rep(1:2, each = n_samples),
      do.call(cbind, variables)
    )
  }
  
  #' Compute U-test p-values for each variable distinguishing classes
  calculate_Utest_pvalues <- function(Data, Cls) {
    feature_cols <- setdiff(names(Data), "Cls")
    p_values <- apply(Data[, feature_cols, drop = FALSE], 2, function(x) {
      wilcox.test(x ~ factor(Cls))$p.value
    })
    data.frame(Variable_No = seq_along(p_values), p_values = p_values)
  }
  
  #' Plot significance comparison across multiple downsampling iterations
  plot_significance_comparison_multi <- function(original_pvals, reduced_pvals_list, removed_pvals_list) {
    original_data <- data.frame(
      Variable_No = original_pvals$Variable_No,
      neg_log_p = -log10(original_pvals$p_values),
      Dataset = "Original",
      Iteration = "Original"
    )
    
    reduced_data <- do.call(rbind, lapply(seq_along(reduced_pvals_list), function(i) {
      data.frame(
        Variable_No = reduced_pvals_list[[i]]$Variable_No,
        neg_log_p = -log10(reduced_pvals_list[[i]]$p_values),
        Dataset = "Reduced",
        Iteration = paste0("Iter_", i)
      )
    }))
    
    removed_data <- do.call(rbind, lapply(seq_along(removed_pvals_list), function(i) {
      data.frame(
        Variable_No = removed_pvals_list[[i]]$Variable_No,
        neg_log_p = -log10(removed_pvals_list[[i]]$p_values),
        Dataset = "Removed",
        Iteration = paste0("Iter_", i)
      )
    }))
    
    combined_data <- rbind(original_data, reduced_data, removed_data)
    combined_data$Dataset <- factor(combined_data$Dataset, levels = c("Original", "Reduced", "Removed"))
    
    plot_theme <- theme_light() +
      theme(
        legend.position = "bottom",
        legend.background = element_rect(fill = alpha("white", 0.5)),
        strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black"),
        panel.grid.minor = element_blank()
      )
    
    ggplot(combined_data, aes(x = Variable_No, y = neg_log_p)) +
      geom_point(data = subset(combined_data, Dataset == "Original"), color = "steelblue", alpha = 0.8, size = 1.5) +
      geom_point(data = subset(combined_data, Dataset == "Reduced"), aes(color = Dataset), alpha = 0.6, size = 1.2) +
      geom_point(data = subset(combined_data, Dataset == "Removed"), aes(color = Dataset), alpha = 0.6, size = 1.2) +
      geom_hline(yintercept = -log10(significance_threshold), color = "red", linetype = "dashed", size = 1) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 3) +
      scale_color_colorblind() +
      plot_theme +
      labs(
        x = "Variable Number",
        y = "-log10(p-value)",
        title = paste0("Statistical Significance Comparison: Original vs ", length(reduced_pvals_list), " Downsampling Iterations"),
        subtitle = paste("Based on", format(nTrials, scientific = FALSE), "trials per iteration;",
                         "CheckRemoved =", CheckRemoved, ", OptimizeBetween =", OptimizeBetween, ", NonNoiseSelection =", NonNoiseSelection),
        color = "Dataset Type"
      )
  }
  
  
  #' ===== Main Experiment Execution =====
  
  set.seed(random_seed)
  
  cat("Generating synthetic dataset with noise variables...\n")
  base_data <- generate_synthetic_data(samples_per_class, num_variables, random_seed)

  # Generate noise variables matching value range of informative features
  noise_data <- lapply(1:num_variables, function(i) {
    class1_mean <- seq(1, 50, length.out = num_variables)[i]
    class2_mean <- seq(1, 48, length.out = num_variables)[i]
    overall_mean <- mean(c(class1_mean, class2_mean))
    sd <- 5

    range_min <- overall_mean - 3 * sd
    range_max <- overall_mean + 3 * sd

    runif(2 * samples_per_class, min = range_min, max = range_max)
  })

  # Generate uniform variables with ascending significant differences between classes
  uniform_class_diff_data <- lapply(1:num_variables, function(i) {
    # Create ascending difference between classes (from very subtle to small)
    # Variable 1 has smallest difference, variable 50 has largest difference
    base_range <- 100  # Base range for uniform distribution

    # Very small ascending difference: from 0.05% to 2% of base range
    class_difference <- seq(0.0005, 0.025, length.out = num_variables)[i] * base_range

    # Set class centers with increasing separation
    class1_center <- 50  # Fixed center for class 1
    class2_center <- class1_center + class_difference  # Increasing separation for class 2

    # Width of uniform distribution for each class (fixed)
    class_width <- base_range * 0.2  # 20% of base range

    # Generate uniform data for each class
    class1_data <- runif(samples_per_class,
                         min = class1_center - class_width/2,
                         max = class1_center + class_width/2)

    class2_data <- runif(samples_per_class,
                         min = class2_center - class_width/2,
                         max = class2_center + class_width/2)

    c(class1_data, class2_data)
  })

  # Combine all data types
  final_data <- cbind(base_data,
                      do.call(cbind, noise_data),
                      do.call(cbind, uniform_class_diff_data))

  # Update column names to reflect the different variable types
  total_vars <- num_variables + num_variables + num_variables
  names(final_data) <- c("Cls",
                         paste0("GMMVar", 1:num_variables),
                         paste0("NoiseVar", 1:num_variables),
                         paste0("UniformVar", 1:num_variables))

  data_df <- list(Data = final_data[, -1], Cls = final_data$Cls)
  
  cat("Calculating p-values for original dataset...\n")
  pval_results <- calculate_Utest_pvalues(final_data, final_data$Cls)
  
  cat(paste0("Performing ", nSamples, " downsampling iterations...\n"))
  downsampled_results <- lapply(1:nSamples, function(i) {
    cat(paste0("  Iteration ", i, "/", nSamples, "\n"))
    result <- opdisDownsampling::opdisDownsampling(
      Data = data_df$Data,
      Cls = data_df$Cls,
      Size = downsampling_size,
      Seed = i + (i - 1) * 1e6,
      nTrials = nTrials,
      CheckRemoved = CheckRemoved,
      OptimizeBetween = OptimizeBetween,
      MaxCores = min(nTrials, parallel::detectCores()-1),
      TestStat = TestStat, 
      PCAimportance = FALSE,
      NonNoiseSelection = NonNoiseSelection
    )
    list(
      reduced_data = cbind(Cls = result$ReducedData$Cls, result$ReducedData[, -ncol(result$ReducedData)]),
      removed_data = cbind(Cls = result$RemovedData$Cls, result$RemovedData[, -ncol(result$RemovedData)])
    )
  })
  
  names(downsampled_results) <- paste0("Iteration",1:nSamples)
  
  cat("Calculating p-values for reduced datasets...\n")
  pval_reduced_list <- lapply(downsampled_results, function(res) {
    calculate_Utest_pvalues(res$reduced_data, res$reduced_data$Cls)
  })
  
  cat("Calculating p-values for removed datasets...\n")
  pval_removed_list <- lapply(downsampled_results, function(res) {
    calculate_Utest_pvalues(res$removed_data, res$removed_data$Cls)
  })
  
  cat("Creating visualization...\n")
  comparison_plot <- plot_significance_comparison_multi(pval_results, pval_reduced_list, pval_removed_list)
  
  # Extract p-values matrices for Kendall's tau
  pval_reduced_matrix <- do.call(cbind, lapply(pval_reduced_list, `[[`, "p_values"))
  pval_removed_matrix <- do.call(cbind, lapply(pval_removed_list, `[[`, "p_values"))
  
  # Compute Kendall's tau correlations
  tau_reduced <- apply(pval_reduced_matrix, 2, cor, y = pval_results$p_values, method = "kendall")
  tau_removed <- apply(pval_removed_matrix, 2, cor, y = pval_results$p_values, method = "kendall")
  tau_reduced_vs_removed <- mapply(cor, as.data.frame(pval_reduced_matrix), as.data.frame(pval_removed_matrix), MoreArgs = list(method = "kendall"))
  
  # Prepare data for plotting correlations
  correlation_data <- rbind(
    data.frame(Iteration = seq_along(tau_reduced), Tau = tau_reduced, Comparison = "Reduced vs Original"),
    data.frame(Iteration = seq_along(tau_removed), Tau = tau_removed, Comparison = "Removed vs Original"),
    data.frame(Iteration = seq_along(tau_reduced_vs_removed), Tau = tau_reduced_vs_removed, Comparison = "Reduced vs Removed")
  )
  
  # Summary statistics for annotations
  median_tau_reduced <- median(tau_reduced)
  median_tau_removed <- median(tau_removed)
  median_tau_reduced_vs_removed <- median(tau_reduced_vs_removed)
  
  variance_tau_reduced <- var(tau_reduced)
  variance_tau_removed <- var(tau_removed)
  variance_tau_reduced_vs_removed <- var(tau_reduced_vs_removed)
  
  min_tau_reduced <- min(tau_reduced)
  min_tau_removed <- min(tau_removed)
  min_tau_reduced_vs_removed <- min(tau_reduced_vs_removed)
  
  # Plot Kendall's tau correlations as dot plots with overlaid box plots and statistics as annotations
  # Prepare annotation data
  annotation_df <- data.frame(
    Comparison = c("Reduced vs Original", "Removed vs Original", "Reduced vs Removed"),
    x = c("Reduced vs Original", "Removed vs Original", "Reduced vs Removed"),
    y = c(1.2, 1.2, 1.2),
    label = c(
      paste("Med:", round(median_tau_reduced, 3), "\nVar:", round(variance_tau_reduced, 3), "\nMin:", round(min_tau_reduced, 3)),
      paste("Med:", round(median_tau_removed, 3), "\nVar:", round(variance_tau_removed, 3), "\nMin:", round(min_tau_removed, 3)),
      paste("Med:", round(median_tau_reduced_vs_removed, 3), "\nVar:", round(variance_tau_reduced_vs_removed, 3), "\nMin:", round(min_tau_reduced_vs_removed, 3))
    )
  )
  
  p_kendall_correlations_box <- ggplot(correlation_data, aes(x = factor(Comparison), y = Tau, fill = Comparison)) +
    geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk") +
    geom_point(position = position_jitter(width = 0.05), alpha = 1) +
    geom_text(
      data = annotation_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1, vjust = .5, size = 3.5, angle = 90, fontface = "plain", color = "darkblue"
    ) +
    theme_light() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Comparison",
      y = "Kendall's Tau",
      title = "Kendall's Tau ",
      subtitle = "p (group differences)"
    )
  
  # Combine and show plots
  combined_plot <- cowplot::plot_grid(comparison_plot, p_kendall_correlations_box,
                                      labels = "AUTO", rel_widths = c(4, 1))
  print(combined_plot)
  
  # Save plot to file
  ggsave(paste0("combined_plot_", 
                format(nTrials, scientific = FALSE), "trials_per_iteration_CheckRemoved", CheckRemoved, "_OptimizeBetween", OptimizeBetween,
                "_NonNoiseSelection", NonNoiseSelection, ".svg"),
         combined_plot, width = 14, height = 6)
  
  # Print summary statistics to console
  print(paste("Reduced vs Original - Med:", round(median_tau_reduced, 3), "Var:", round(variance_tau_reduced, 3), "Min:", round(min_tau_reduced, 3)))
  print(paste("Removed vs Original - Med:", round(median_tau_removed, 3), "Var:", round(variance_tau_removed, 3), "Min:", round(min_tau_removed, 3)))
  print(paste("Reduced vs Removed - Med:", round(median_tau_reduced_vs_removed, 3), "Var:", round(variance_tau_reduced_vs_removed, 3), "Min:", round(min_tau_reduced_vs_removed, 3)))
  
  cat("\nAnalysis complete for current parameter set.\n\n")
  
  return(list(p_values_orig = pval_results,
              pval_reduced_list = pval_reduced_list,
              pval_removed_list = pval_removed_list,
              downsampled_data = downsampled_results,
              correlations = correlation_data,
              combined_plot = combined_plot, 
              comparison_plot = comparison_plot, 
              p_kendall_correlations_box = p_kendall_correlations_box
  ))
}

# Generate all parameter combinations for data splitting options
params <- expand.grid(
  CheckRemoved = c(TRUE, FALSE),
  OptimizeBetween = c(TRUE, FALSE),
  NonNoiseSelection = c(TRUE, FALSE)
)


# Conditionally run experiment either once or for all parameter combinations
if (nTrials == 1) {
  # Single run with default parameters (all FALSE)
  results_experiment_ascending_significance_1Trial <- run_experiment(FALSE, FALSE, FALSE, FALSE)
} else if (nTrials > 1) {
  # Run all parameter combinations
  results_experiment_ascending_significance_nTrials <-
    lapply(seq_len(nrow(params)), function(i){
      res <- run_experiment(
        params$CheckRemoved[i],
        params$OptimizeBetween[i],
        params$NonNoiseSelection[i]
      )
      return(res)
    })

  # Name the results list with parameter combinations for easy identification
  names(results_experiment_ascending_significance_nTrials) <- paste0(
    "CR", params$CheckRemoved,
    "_OB", params$OptimizeBetween,
    "_NNS", params$NonNoiseSelection
  )
}

cat("All analyses complete!\n")


# Extract correlations from single trial experiment
single_trial_correlations <- results_experiment_ascending_significance_1Trial$correlations

# Add parameter combination identifier (since it's the default: all FALSE)
single_trial_correlations$ParameterCombination <- "CRFALSE_CTFALSE_OBFALSE_NNSFALSE"

# Extract all correlations data frames and combine them
all_correlations <- do.call(rbind, lapply(names(results_experiment_ascending_significance_nTrials), function(param_name) {
  correlations_df <- results_experiment_ascending_significance_nTrials[[param_name]]$correlations
  # Add parameter combination as a column
  correlations_df$ParameterCombination <- param_name
  return(correlations_df)
}))

# Reset row names to avoid duplicates
rownames(all_correlations) <- NULL

# View the structure of the combined data frame
str(all_correlations)
head(all_correlations)

# Parse parameter names into separate columns
library(tidyr)
all_correlations <- all_correlations %>%
  separate(ParameterCombination,
           into = c("CheckRemoved", "OptimizeBetween", "NonNoiseSelection"),
           sep = "_",
           remove = FALSE) %>%
  mutate(
    CheckRemoved = gsub("CR", "", CheckRemoved),
    OptimizeBetween = gsub("OB", "", OptimizeBetween),
    NonNoiseSelection = gsub("NNS", "", NonNoiseSelection)
  )


# Create the styled plot for single trial
p_single_trial_correlations <- ggplot(data = single_trial_correlations, aes(x = Comparison, y = Tau)) +
  geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk") +
  geom_point(aes(color = as.factor(Iteration)), position = position_jitter(width = 0.05), alpha = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "salmon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3") +
  scale_color_colorblind(name = "Iteration") +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0.5)),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(colour = "black"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(color = guide_legend(nrow = 1)) +
  labs(
    x = "Comparison",
    y = "Kendall's Tau",
    title = "Single Trial Experiment",
    subtitle = paste("1 trial per iteration,", nSamples, "iterations"),
    color = "Iteration"
  )

# Create the styled plot for multi-trial experiment
p_all_correlations <- ggplot(data = all_correlations, aes(x = Comparison, y = Tau)) +
  geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk") +
  geom_point(aes(color = as.factor(Iteration)), position = position_jitter(width = 0.05), alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "salmon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3") +
  facet_wrap(.~ParameterCombination) +
  scale_color_colorblind(name = "Iteration") +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0.5)),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(colour = "black"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(color = guide_legend(nrow = 1)) +
  labs(
    x = "Comparison",
    y = "Kendall's Tau",
    title = "Multiple Parameter Combinations",
    subtitle = paste(format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations"),
    color = "Iteration"
  )

# Combine the two plots side by side
p_combined_correlations <- cowplot::plot_grid(
  p_single_trial_correlations,
  p_all_correlations,
  labels = c("A", "B"),
  rel_widths = c(1, 4),
  align = "h",  axis = "b",
  ncol = 2
)

# Display the combined plot
print(p_combined_correlations)

# Save the combined plot
ggsave(paste0("combined_correlations_comparison_",
              format(nTrials, scientific = FALSE), "trials_",
              nSamples, "iterations.svg"),
       p_combined_correlations, width = 15, height = 15)


# Check for identical parameter combinations in all_correlations
check_identical_parameter_combinations <- function(df) {
  # Get unique parameter combinations
  param_combinations <- unique(df$ParameterCombination)

  # Extract Tau values for each parameter combination
  tau_blocks <- lapply(param_combinations, function(param) {
    subset_df <- df[df$ParameterCombination == param, ]
    # Sort by iteration to ensure correct order
    subset_df <- subset_df[order(subset_df$Iteration), ]
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

    for (j in (i+1):length(tau_blocks)) {
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
    # Extract individual parameters from string like "CRFALSE_CTFALSE_OBFALSE_NNSFALSE"
    parts <- strsplit(param_combo, "_")[[1]]
    params <- list()

    for (part in parts) {
      if (grepl("^CR", part)) {
        params$CheckRemoved <- grepl("TRUE", part)
      } else if (grepl("^OB", part)) {
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
    param_names <- c("CheckRemoved", "OptimizeBetween", "NonNoiseSelection")
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

# Additional function for custom single runs
run_custom_experiment <- function(CheckRemoved = FALSE, 
                                  OptimizeBetween = FALSE, NonNoiseSelection = FALSE,
                                  custom_nTrials = NULL) {
  
  if (is.null(custom_nTrials)) {
    custom_nTrials <- nTrials  # Use global nTrials if not specified
  }
  
  cat(paste0("Running custom experiment with nTrials = ", custom_nTrials,
             " and parameters: CR=", CheckRemoved, 
             ", OB=", OptimizeBetween, ", NNS=", NonNoiseSelection, "\n"))
  
  # Temporarily store original nTrials and set custom value
  original_nTrials <- nTrials
  assign("nTrials", custom_nTrials, envir = .GlobalEnv)
  
  # Run the experiment
  result <- run_experiment(CheckRemoved, OptimizeBetween, NonNoiseSelection)
  
  # Restore original nTrials
  assign("nTrials", original_nTrials, envir = .GlobalEnv)
  
  return(result)
}

# run_custom_experiment(custom_nTrials = 10, NonNoiseSelection = T)