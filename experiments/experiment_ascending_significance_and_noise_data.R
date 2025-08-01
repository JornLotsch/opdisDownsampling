# ===============================================================================
# ARTIFICIAL DATA SETS GENERATION AND EVALUATION SCRIPT
# ===============================================================================
#
# Purpose: Generates synthetic datasets for evaluating data splitting methods,
#          featuring variables with clear class separation and added noise variables.
#          Performs assessment of data splitting parameter effects on feature
#          significance preservation.
#
# This is a demonstration script - works as-is but not optimized for production use.
# ===============================================================================

# Attempt to set working directory to script location (for RStudio)
tryCatch({
  if (exists("rstudioapi::getSourceEditorContext", mode = "function")) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}, error = function(e) {
  message("Unable to set working directory automatically. Please set it manually if needed.")
})

setwd("/home/joern/Aktuell/DownSamplingStructure/12RLibrary/opdisDownsampling/experiments")

# ===============================================================================
# LOAD REQUIRED LIBRARIES
# ===============================================================================
library(stats)      # Statistical functions
library(ggplot2)    # Data visualization
library(ggthemes)   # Additional ggplot2 themes
library(dplyr)      # Data manipulation
library(parallel)      # Parallel processing
library(pbmcapply)      # Progress bars for parallel processing
library(cowplot)      # Plot arrangement
library(reshape2)      # Data reshaping
library(tidyr)      # Data tidying

# ===============================================================================
# EXPERIMENTAL CONFIGURATION PARAMETERS
# ===============================================================================
random_seed <- 42            # Random seed for reproducibility
samples_per_class <- 100     # Number of samples per class
num_variables <- 50          # Number of variables (features) per data type
significance_threshold <- 0.05  # Statistical significance threshold for p-values
downsampling_size <- 0.8     # Proportion of original data to retain (80%)
nSamples <- 10               # Number of downsampling iterations to perform
nTrials <- 100000                # Number of optimization trials per downsampling iteration
TestStat = "ad"              # Statistical test for distribution comparison (Anderson-Darling)
use_ABC_for_Pvalues = TRUE   # Calculate Tau only for relevant features
WorstSample = FALSE          # Test if results worsen when reversing the selection criterions

# Default values for data splitting optimization options
OptimizeBetween <- FALSE     # Whether to optimize between reduced and removed sets
NonNoiseSelection <- FALSE   # Whether to use noise detection for variable pre-selection


# ===============================================================================
# MAIN EXPERIMENT FUNCTION
# ===============================================================================
# Function to run experiment given parameters
run_experiment <- function(OptimizeBetween, NonNoiseSelection) {

  # Print current parameter combination for tracking progress
  print(c(OptimizeBetween, NonNoiseSelection))

  # ---------------------------------------------------------------------------
  # DATA GENERATION FUNCTIONS
  # ---------------------------------------------------------------------------

    #' Generate synthetic dataset with ascending class separation
    #' Creates Gaussian-distributed variables where class separation increases progressively
  generate_synthetic_data <- function(n_samples, n_vars, seed) {
    # Define class means with increasing separation between classes
    class1_means <- seq(1, 50, length.out = n_vars)
    class2_means <- seq(1, 48, length.out = n_vars)

    # Generate variables for both classes with normal distributions
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
    #' Uses Wilcoxon rank-sum test to assess class discrimination ability
  calculate_Utest_pvalues <- function(Data, Cls) {
    feature_cols <- setdiff(names(Data), "Cls")
    p_values <- apply(Data[, feature_cols, drop = FALSE], 2, function(x) {
      wilcox.test(x ~ factor(Cls))$p.value
    })
    data.frame(Variable_No = seq_along(p_values), p_values = p_values)
  }

    #' Plot significance comparison across multiple downsampling iterations
    #' Creates faceted visualization comparing original, reduced, and removed datasets
  plot_significance_comparison_multi <- function(original_pvals, reduced_pvals_list, removed_pvals_list) {
    # Prepare original dataset p-values for plotting
    original_data <- data.frame(
      Variable_No = original_pvals$Variable_No,
      neg_log_p = -log10(original_pvals$p_values),
      Dataset = "Original",
      Iteration = "Original"
    )

    # Prepare reduced dataset p-values from all iterations
    reduced_data <- do.call(rbind, lapply(seq_along(reduced_pvals_list), function(i) {
      data.frame(
        Variable_No = reduced_pvals_list[[i]]$Variable_No,
        neg_log_p = -log10(reduced_pvals_list[[i]]$p_values),
        Dataset = "Reduced",
        Iteration = paste0("Iter_", i)
      )
    }))

    # Prepare removed dataset p-values from all iterations
    removed_data <- do.call(rbind, lapply(seq_along(removed_pvals_list), function(i) {
      data.frame(
        Variable_No = removed_pvals_list[[i]]$Variable_No,
        neg_log_p = -log10(removed_pvals_list[[i]]$p_values),
        Dataset = "Removed",
        Iteration = paste0("Iter_", i)
      )
    }))

    # Combine all data and set factor levels for consistent plotting
    combined_data <- rbind(original_data, reduced_data, removed_data)
    combined_data$Dataset <- factor(combined_data$Dataset, levels = c("Original", "Reduced", "Removed"))

    # Define plot theme for consistent styling
    plot_theme <- theme_light() +
      theme(
        legend.position = "bottom",
        legend.background = element_rect(fill = alpha("white", 0.5)),
        strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black"),
        panel.grid.minor = element_blank()
      )

    # Create the multi-panel plot
    ggplot(combined_data, aes(x = Variable_No, y = neg_log_p, shape = as.factor(Iteration))) +
      geom_point(data = subset(combined_data, Dataset == "Original"), color = "steelblue", alpha = 0.8, size = 1.5) +
      geom_point(data = subset(combined_data, Dataset == "Reduced"), aes(color = Dataset), alpha = 0.6, size = 1.2) +
      geom_point(data = subset(combined_data, Dataset == "Removed"), aes(color = Dataset), alpha = 0.6, size = 1.2) +
      geom_hline(yintercept = -log10(significance_threshold), color = "salmon", linetype = "dashed", size = 1) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 3) +
      scale_color_colorblind() +
      plot_theme +
      labs(
        x = "Variable Number",
        y = "-log10(p-value)",
        title = paste0("Statistical Significance Comparison: Original vs ", length(reduced_pvals_list), " Downsampling Iterations"),
        subtitle = paste0("Based on", format(nTrials, scientific = FALSE), "trials per iteration; ","OptimizeBetween =", OptimizeBetween, ", NonNoiseSelection =", NonNoiseSelection),
        color = "Dataset Type"
      ) +
      guides(shape = "none")
  }

  # ---------------------------------------------------------------------------
  # MAIN EXPERIMENT EXECUTION
  # ---------------------------------------------------------------------------

  set.seed(random_seed)

  cat("Generating synthetic dataset with noise variables...\n")
  # Generate base Gaussian mixture data with ascending class separation
  base_data <- generate_synthetic_data(samples_per_class, num_variables, random_seed)

  # Generate pure noise variables (uniform random within same value range as informative features)
  noise_data <- lapply(1:num_variables, function(i) {
    class1_mean <- seq(1, 50, length.out = num_variables)[i]
    class2_mean <- seq(1, 48, length.out = num_variables)[i]
    overall_mean <- mean(c(class1_mean, class2_mean))
    sd <- 5

    # Calculate range based on 3 standard deviations
    range_min <- overall_mean - 3 * sd
    range_max <- overall_mean + 3 * sd

    # Generate uniform random values within this range
    runif(2 * samples_per_class, min = range_min, max = range_max)
  })

  # Generate uniform variables with very subtle but ascending significant differences between classes
  uniform_class_diff_data <- lapply(1:num_variables, function(i) {
    # Create ascending difference between classes (from very subtle to small)
    # Variable 1 has smallest difference, variable 50 has largest difference
    base_range <- 100  # Base range for uniform distribution

    # Very small ascending difference: from 0.05% to 2% of base range
    class_difference <- seq(0.0005, 0.025, length.out = num_variables)[i] * base_range

    # Set class centers with increasing separation
    class1_center <- 50  # Fixed center for class 1
    class2_center <- class1_center + class_difference  # Increasing separation for class 2

    # Width of uniform distribution for each class (fixed at 20% of base range)
    class_width <- base_range * 0.2

    # Generate uniform data for each class
    class1_data <- runif(samples_per_class,
                         min = class1_center - class_width/2,
                         max = class1_center + class_width/2)

    class2_data <- runif(samples_per_class,
                         min = class2_center - class_width/2,
                         max = class2_center + class_width/2)

    c(class1_data, class2_data)
  })

  # Combine all three data types into final dataset
  final_data <- cbind(base_data,
                      do.call(cbind, noise_data),
                      do.call(cbind, uniform_class_diff_data))

  # Update column names to reflect the different variable types for clarity
  names(final_data) <- c("Cls",
                         paste0("GMMVar", 1:num_variables),      # Gaussian Mixture Model variables
                         paste0("NoiseVar", 1:num_variables),    # Pure noise variables
                         paste0("UniformVar", 1:num_variables))  # Uniform variables with class differences

  # Prepare data structure for opdisDownsampling function
  data_df <- list(Data = final_data[, -1], Cls = final_data$Cls)

  cat("Calculating p-values for original dataset...\n")
  # Calculate statistical significance for all variables in original dataset
  pval_results <- calculate_Utest_pvalues(final_data, final_data$Cls)

  cat(paste0("Performing ", nSamples, " downsampling iterations...\n"))
  # Perform multiple downsampling iterations with different seeds
  downsampled_results <- lapply(1:nSamples, function(i) {
    cat(paste0("  Iteration ", i, "/", nSamples, "\n"))
    result <- opdisDownsampling::opdisDownsampling(
      Data = data_df$Data,
      Cls = data_df$Cls,
      Size = downsampling_size,
      Seed = i + (i - 1) * 1e6,       
      nTrials = nTrials,
      OptimizeBetween = OptimizeBetween,    
      MaxCores = min(nTrials, parallel::detectCores()-1),
      TestStat = TestStat,
      PCAimportance = FALSE,
      NonNoiseSelection = NonNoiseSelection,
      WorstSample = WorstSample              
    )
    # Return both reduced (kept) and removed data with class labels
    list(
      reduced_data = cbind(Cls = result$ReducedData$Cls, result$ReducedData[, -ncol(result$ReducedData)]),
      removed_data = cbind(Cls = result$RemovedData$Cls, result$RemovedData[, -ncol(result$RemovedData)])
    )
  })

  names(downsampled_results) <- paste0("Iteration",1:nSamples)

  cat("Calculating p-values for reduced datasets...\n")
  # Calculate p-values for each reduced dataset
  pval_reduced_list <- lapply(downsampled_results, function(res) {
    calculate_Utest_pvalues(res$reduced_data, res$reduced_data$Cls)
  })

  cat("Calculating p-values for removed datasets...\n")
  # Calculate p-values for each removed dataset
  pval_removed_list <- lapply(downsampled_results, function(res) {
    calculate_Utest_pvalues(res$removed_data, res$removed_data$Cls)
  })

  cat("Creating visualization...\n")
  # Create main significance comparison visualization
  comparison_plot <- plot_significance_comparison_multi(pval_results, pval_reduced_list, pval_removed_list)

  # ---------------------------------------------------------------------------
  # JACCARD INDEX ANALYSIS (Feature Selection Overlap)
  # ---------------------------------------------------------------------------

  # Calculate Jaccard Index for feature selection overlap between sets
  calculate_jaccard_index <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    if (union == 0) return(1) # Both sets are empty - perfect agreement
    return(intersection / union)
  }

  # Enhanced p-value calculation with three-fold Jaccard analysis
  # Compares feature selection overlap between: Original-Reduced, Original-Removed, Reduced-Removed
  calculate_pvals_with_jaccard_threefold <- function(pval_results, pval_reduced_list, pval_removed_list,
                                                     significance_threshold = 0.05) {

    # Apply FDR correction instead of Bonferroni for multiple testing
    n_variables <- nrow(pval_results)
    fdr_adjusted_pvals <- p.adjust(pval_results$p_values, method = "fdr")

    cat("Calculating three-fold Jaccard indices for feature selection overlap...\n")

    # Identify significant features for original data using both raw and FDR-corrected p-values
    original_significant_raw <- which(pval_results$p_values < significance_threshold)
    original_significant_fdr <- which(fdr_adjusted_pvals < significance_threshold)

    # Initialize results data frame for all three comparison types
    jaccard_results <- data.frame(
      Iteration = integer(),
      Comparison = character(),
      Jaccard_Raw = numeric(),
      Jaccard_FDR = numeric(),
      stringsAsFactors = FALSE
    )

    # Process each downsampling iteration
    for (i in 1:length(pval_reduced_list)) {
      # Get significant features for this iteration (both reduced and removed sets)
      reduced_significant_raw <- which(pval_reduced_list[[i]]$p_values < significance_threshold)
      reduced_fdr_adjusted <- p.adjust(pval_reduced_list[[i]]$p_values, method = "fdr")
      reduced_significant_fdr <- which(reduced_fdr_adjusted < significance_threshold)

      removed_significant_raw <- which(pval_removed_list[[i]]$p_values < significance_threshold)
      removed_fdr_adjusted <- p.adjust(pval_removed_list[[i]]$p_values, method = "fdr")
      removed_significant_fdr <- which(removed_fdr_adjusted < significance_threshold)

      # 1. Reduced vs Original - How well does downsampling preserve significant features?
      jaccard_raw_red_orig <- calculate_jaccard_index(original_significant_raw, reduced_significant_raw)
      jaccard_fdr_red_orig <- calculate_jaccard_index(original_significant_fdr, reduced_significant_fdr)

      jaccard_results <- rbind(jaccard_results, data.frame(
        Iteration = i,
        Comparison = "Reduced vs Original",
        Jaccard_Raw = jaccard_raw_red_orig,
        Jaccard_FDR = jaccard_fdr_red_orig,
        stringsAsFactors = FALSE
      ))

      # 2. Removed vs Original - Do removed samples contain different significant features?
      jaccard_raw_rem_orig <- calculate_jaccard_index(original_significant_raw, removed_significant_raw)
      jaccard_fdr_rem_orig <- calculate_jaccard_index(original_significant_fdr, removed_significant_fdr)

      jaccard_results <- rbind(jaccard_results, data.frame(
        Iteration = i,
        Comparison = "Removed vs Original",
        Jaccard_Raw = jaccard_raw_rem_orig,
        Jaccard_FDR = jaccard_fdr_rem_orig,
        stringsAsFactors = FALSE
      ))

      # 3. Reduced vs Removed - How complementary are the feature sets in the two partitions?
      jaccard_raw_red_rem <- calculate_jaccard_index(reduced_significant_raw, removed_significant_raw)
      jaccard_fdr_red_rem <- calculate_jaccard_index(reduced_significant_fdr, removed_significant_fdr)

      jaccard_results <- rbind(jaccard_results, data.frame(
        Iteration = i,
        Comparison = "Reduced vs Removed",
        Jaccard_Raw = jaccard_raw_red_rem,
        Jaccard_FDR = jaccard_fdr_red_rem,
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
  create_jaccard_panel <- function(jaccard_analysis, nSamples, nTrials,
                                   OptimizeBetween, NonNoiseSelection) {

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
    p_jaccard <- ggplot(jaccard_long, aes(x = Comparison, y = Jaccard_Index, color = Comparison)) +
      geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk", outlier.shape = NA, size = 0.3) +
      geom_point(aes(color = Comparison, shape = as.factor(Iteration)), position = position_jitter(width = 0.1), alpha = 0.7, size = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3", size = 0.5) +      # Perfect overlap
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", size = 0.3) +         # Moderate overlap
      geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 0.3, alpha = 0.5) + # No overlap
      facet_wrap(~Correction_Type, scales = "free_x", ncol = 1) +
      scale_color_colorblind() +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
      plot_theme +
      labs(
        x = "Comparison",
        y = "Jaccard Index",
        title = "Feature Selection Overlap"
      ) +
      guides(shape = "none")

    return(p_jaccard)
  }

  # Execute Jaccard analysis workflow
  names(downsampled_results) <- paste0("Iteration", 1:nSamples)

  cat("Calculating p-values for reduced datasets...\n")
  pval_reduced_list <- lapply(downsampled_results, function(res) {
    calculate_Utest_pvalues(res$reduced_data, res$reduced_data$Cls)
  })

  cat("Calculating p-values for removed datasets...\n")
  pval_removed_list <- lapply(downsampled_results, function(res) {
    calculate_Utest_pvalues(res$removed_data, res$removed_data$Cls)
  })

  cat("Calculating three-fold Jaccard indices for feature selection overlap...\n")
  # Perform comprehensive Jaccard analysis
  jaccard_analysis <- calculate_pvals_with_jaccard_threefold(pval_results, pval_reduced_list, pval_removed_list,
                                                             significance_threshold)

  cat("Creating visualizations...\n")
  # Create main significance comparison plot
  comparison_plot <- plot_significance_comparison_multi(pval_results, pval_reduced_list, pval_removed_list)

  # Create the Jaccard index panel
  p_jaccard <- create_jaccard_panel(jaccard_analysis, nSamples, nTrials,
                                    OptimizeBetween, NonNoiseSelection)

  # ---------------------------------------------------------------------------
  # RESULTS SUMMARY AND INTERPRETATION
  # ---------------------------------------------------------------------------

  # Print enhanced summary statistics for Jaccard indices
  cat("\nThree-fold Jaccard Index Summary:\n")
  cat("=================================\n")

  # Calculate summary statistics for each comparison type
  summary_stats <- aggregate(cbind(Jaccard_Raw, Jaccard_FDR) ~ Comparison,
                             data = jaccard_analysis$jaccard_results,
                             FUN = function(x) c(Mean = round(mean(x), 3),
                                                 Median = round(median(x), 3),
                                                 SD = round(sd(x), 3),
                                                 Min = round(min(x), 3),
                                                 Max = round(max(x), 3)))

  print(summary_stats)

  # Additional insights about feature selection
  cat("\nFeature Selection Insights:\n")
  cat("===========================\n")
  cat(paste("Original significant features (raw p < ", significance_threshold, "): ",
            length(jaccard_analysis$original_significant_raw), " out of ", nrow(pval_results), "\n", sep = ""))
  cat(paste("Original significant features (FDR corrected): ",
            length(jaccard_analysis$original_significant_fdr), " out of ", nrow(pval_results), "\n", sep = ""))

  # Interpretation guidelines for users
  cat("\nInterpretation Guidelines:\n")
  cat("=========================\n")
  cat("• Reduced vs Original: How well does downsampling preserve significant features?\n")
  cat("• Removed vs Original: Do removed samples contain different significant features?\n")
  cat("• Reduced vs Removed: How complementary are the feature sets in the two partitions?\n")
  cat("• High Jaccard (>0.7): Strong overlap\n")
  cat("• Medium Jaccard (0.3-0.7): Moderate overlap\n")
  cat("• Low Jaccard (<0.3): Weak overlap\n")

  # ---------------------------------------------------------------------------
  # CORRELATION ANALYSIS (Kendall's Tau)
  # ---------------------------------------------------------------------------

  # Enhanced validation and error handling
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
  if (use_ABC_for_Pvalues) {
    
    tryCatch({
      feature_indices_orig_AB <- ABCanalysis::ABCanalysis(-log10(pval_orig))$Cind
    }, error = function(e) {
      stop(paste("ABCanalysis failed:", e$message))
    })

    if (length(feature_indices_orig_AB) == 0) {
      warning("No features selected by ABCanalysis - skipping filtering")
    } else {
      if (any(feature_indices_orig_AB > nrow(pval_orig)) || any(feature_indices_orig_AB < 1)) {
        stop("Invalid feature indices detected")
      }

      # Check remaining dimensions after filtering
      remaining_features <- length(pval_orig) - length(feature_indices_orig_AB)
      if (remaining_features <= 1) {
        stop("Too few features remaining after ABC filtering")
      }

      # Keep matrix structure for correlation analysis
      if (is.matrix(pval_orig)) {
        pval_orig <- pval_orig[-feature_indices_orig_AB, , drop = FALSE]
      } else {
        pval_orig <- pval_orig[-feature_indices_orig_AB]
      }
      pval_removed_matrix <- pval_removed_matrix[-feature_indices_orig_AB, , drop = FALSE]
      pval_reduced_matrix <- pval_reduced_matrix[-feature_indices_orig_AB, , drop = FALSE]
    }
  }

  # Ensure pval_orig is a vector for correlation analysis
  pval_orig <- as.vector(pval_orig)

  # Compute Kendall's tau correlations with error handling
  tau_reduced <- tryCatch({
    apply(pval_reduced_matrix, 2, function(x) {
      cor(x, pval_orig, method = "kendall", use = "complete.obs")
    })
  }, error = function(e) {
    warning(paste("Error computing tau_reduced:", e$message))
    rep(NA, ncol(pval_reduced_matrix))
  })

  tau_removed <- tryCatch({
    apply(pval_removed_matrix, 2, function(x) {
      cor(x, pval_orig, method = "kendall", use = "complete.obs")
    })
  }, error = function(e) {
    warning(paste("Error computing tau_removed:", e$message))
    rep(NA, ncol(pval_removed_matrix))
  })

  tau_reduced_vs_removed <- tryCatch({
    mapply(function(x, y) {
      cor(x, y, method = "kendall", use = "complete.obs")
    }, as.data.frame(pval_reduced_matrix), as.data.frame(pval_removed_matrix))
  }, error = function(e) {
    warning(paste("Error computing tau_reduced_vs_removed:", e$message))
    rep(NA, ncol(pval_reduced_matrix))
  })

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
  
  # Calculate summary statistics for annotations
  median_tau_reduced <- median(tau_reduced)
  median_tau_removed <- median(tau_removed)
  median_tau_reduced_vs_removed <- median(tau_reduced_vs_removed)

  variance_tau_reduced <- var(tau_reduced)
  variance_tau_removed <- var(tau_removed)
  variance_tau_reduced_vs_removed <- var(tau_reduced_vs_removed)

  min_tau_reduced <- min(tau_reduced)
  min_tau_removed <- min(tau_removed)
  min_tau_reduced_vs_removed <- min(tau_reduced_vs_removed)

  # Create annotation data frame for plot statistics
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

  # Create Kendall's tau correlation plot with box plots and individual points
  p_kendall_correlations_box <- ggplot(correlation_data, aes(x = factor(Comparison), y = Tau, fill = Comparison)) +
    geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk") +
    geom_point(aes(shape = as.factor(Iteration)), position = position_jitter(width = 0.05), alpha = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon") +     # No correlation line
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3") + # Perfect correlation line
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
    ) +
    guides(shape = "none")

  # ---------------------------------------------------------------------------
  # FINAL VISUALIZATION AND OUTPUT
  # ---------------------------------------------------------------------------

  # Combine and display all plots in a comprehensive layout
  combined_plot <- cowplot::plot_grid(comparison_plot, p_kendall_correlations_box, p_jaccard,
                                      labels = "AUTO", rel_widths = c(4, 1, 1), nrow = 1)
  print(combined_plot)

  # Save combined plot to file with descriptive filename
  ggsave(paste0("combined_plot_",
                format(nTrials, scientific = FALSE), "trials_per_iteration", "_OptimizeBetween", OptimizeBetween,
                "_NonNoiseSelection", NonNoiseSelection, ".svg"),
         combined_plot, width = 16, height = 12)

  # Print summary statistics to console for quick reference
  print(paste("Reduced vs Original - Med:", round(median_tau_reduced, 3), "Var:", round(variance_tau_reduced, 3), "Min:", round(min_tau_reduced, 3)))
  print(paste("Removed vs Original - Med:", round(median_tau_removed, 3), "Var:", round(variance_tau_removed, 3), "Min:", round(min_tau_removed, 3)))
  print(paste("Reduced vs Removed - Med:", round(median_tau_reduced_vs_removed, 3), "Var:", round(variance_tau_reduced_vs_removed, 3), "Min:", round(min_tau_reduced_vs_removed, 3)))

  cat("\nAnalysis complete for current parameter set.\n\n")

  # Return comprehensive results for further analysis
  return(list(p_values_orig = pval_results,
              pval_reduced_list = pval_reduced_list,
              pval_removed_list = pval_removed_list,
              downsampled_data = downsampled_results,
              correlations = correlation_data,
              combined_plot = combined_plot,
              comparison_plot = comparison_plot,
              p_kendall_correlations_box = p_kendall_correlations_box,
              jaccard_indices = jaccard_analysis$jaccard_results
  ))
}

# ===============================================================================
# CUSTOM SINGLE EXPERIMENT FUNCTION
# ===============================================================================

# Additional function for custom single runs
run_custom_experiment <- function(OptimizeBetween = FALSE, NonNoiseSelection = FALSE,
                                  custom_nTrials = NULL) {

  if (is.null(custom_nTrials)) {
    custom_nTrials <- nTrials  # Use global nTrials if not specified
  }

  cat(paste0("Running custom experiment with nTrials = ", custom_nTrials,
             " and parameters: OB=", OptimizeBetween, ", NNS=", NonNoiseSelection, "\n"))

  # Temporarily store original nTrials and set custom value
  original_nTrials <- nTrials
  assign("nTrials", custom_nTrials, envir = .GlobalEnv)

  # Run the experiment
  result <- run_experiment(OptimizeBetween, NonNoiseSelection)

  # Restore original nTrials
  assign("nTrials", original_nTrials, envir = .GlobalEnv)

  return(result)
}

# Usage example (commented out):
# run_custom_experiment(custom_nTrials = 10, NonNoiseSelection = TRUE)

# ===============================================================================
# EXPERIMENT EXECUTION LOGIC
# ===============================================================================

# Generate all parameter combinations for comprehensive testing
params <- expand.grid(
  OptimizeBetween = c(TRUE, FALSE),      # Test both optimization strategies
  NonNoiseSelection = c(TRUE, FALSE)     # Test both variable selection approaches
)

# Input validation for number of trials
if (!is.numeric(nTrials) || nTrials <= 0) {
  stop("nTrials must be a positive integer")
}

# Execute experiments based on trial configuration
if (nTrials == 1) {
  # Single run with default parameters (all FALSE) - quick test mode
  results_experiment_ascending_significance_1Trial <- run_experiment(FALSE, FALSE)
} else if (nTrials > 1) {
  # Run all parameter combinations - comprehensive analysis mode
  results_experiment_ascending_significance_nTrials <-
    lapply(seq_len(nrow(params)), function(i){
      res <- run_experiment(
        params$OptimizeBetween[i],
        params$NonNoiseSelection[i]
      )
      return(res)
    })

  # Name the results list with parameter combinations for easy identification
  names(results_experiment_ascending_significance_nTrials) <- paste0(
    "_OB", params$OptimizeBetween,
    "_NNS", params$NonNoiseSelection
  )
}

cat("All analyses complete!\n")

# ===============================================================================
# POST-PROCESSING AND COMPARATIVE ANALYSIS
# ===============================================================================

# Extract correlations from single trial experiment if it exists
if (exists("results_experiment_ascending_significance_1Trial")) {
  single_trial_correlations <- results_experiment_ascending_significance_1Trial$correlations
  # Add parameter combination identifier (since it's the default: all FALSE)
  single_trial_correlations$ParameterCombination <- "OBFALSE_NNSFALSE"
} else {
  single_trial_correlations <- data.frame()  # Empty data frame as fallback
}

# Extract all correlations data frames and combine them for comparative analysis
if (exists("results_experiment_ascending_significance_nTrials")) {
  all_correlations <- do.call(rbind, lapply(names(results_experiment_ascending_significance_nTrials), function(param_name) {
    correlations_df <- results_experiment_ascending_significance_nTrials[[param_name]]$correlations
    # Add parameter combination as a column for grouping
    correlations_df$ParameterCombination <- param_name
    return(correlations_df)
  }))

  # Reset row names to avoid duplicates
  rownames(all_correlations) <- NULL

  # Parse parameter names into separate columns for easier analysis
  library(tidyr)
  all_correlations <- all_correlations %>%
    separate(ParameterCombination,
             into = c("OptimizeBetween", "NonNoiseSelection"),
             sep = "_",
             remove = FALSE) %>%
    mutate(
      OptimizeBetween = gsub("OB", "", OptimizeBetween),
      NonNoiseSelection = gsub("NNS", "", NonNoiseSelection)
    )

  # ===============================================================================
  # COMPARATIVE VISUALIZATION
  # ===============================================================================

  # Create comparative plots to visualize differences between parameter settings

  # Single trial correlation plot (if available)
  if (nrow(single_trial_correlations) > 0) {
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
  } else {
    # Create empty plot if no single trial data
    p_single_trial_correlations <- ggplot() + theme_void()
  }

  # Multi-trial correlation plot comparing all parameter combinations
  p_all_correlations <- ggplot(data = all_correlations, aes(x = Comparison, y = Tau)) +
    geom_boxplot(position = "dodge", alpha = 0.5, fill = "cornsilk") +
    geom_point(aes(color = as.factor(Iteration)), position = position_jitter(width = 0.05), alpha = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3") +
    facet_wrap(.~ParameterCombination) +
    scale_color_gdocs(name = "Iteration") +
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

  # Extract and combine Jaccard indices from all experiments for comparative analysis
  all_jaccard_data <- do.call(rbind, lapply(names(results_experiment_ascending_significance_nTrials), function(param_name) {
    jaccard_df <- results_experiment_ascending_significance_nTrials[[param_name]]$jaccard_indices
    # Add parameter combination as a column for grouping
    jaccard_df$ParameterCombination <- param_name
    return(jaccard_df)
  }))

  # Reset row names and parse parameter combinations
  rownames(all_jaccard_data) <- NULL
  all_jaccard_data <- all_jaccard_data %>%
    separate(ParameterCombination,
             into = c("OptimizeBetween", "NonNoiseSelection"),
             sep = "_",
             remove = FALSE) %>%
    mutate(
      OptimizeBetween = gsub("OB", "", OptimizeBetween),
      NonNoiseSelection = gsub("NNS", "", NonNoiseSelection)
    )

  # Create Jaccard comparative visualization
  jaccard_long <- reshape2::melt(all_jaccard_data,
                                 id.vars = c("Iteration", "Comparison", "ParameterCombination", "OptimizeBetween", "NonNoiseSelection"),
                                 measure.vars = c("Jaccard_Raw", "Jaccard_FDR"),
                                 variable.name = "Correction_Type",
                                 value.name = "Jaccard_Index")

  # Clean up labels for better readability
  jaccard_long$Correction_Type <- ifelse(jaccard_long$Correction_Type == "Jaccard_Raw",
                                         "Raw p-values", "FDR corrected")

  # Ensure proper factor ordering for consistent colors across plots
  jaccard_long$Comparison <- factor(jaccard_long$Comparison,
                                    levels = c("Reduced vs Original", "Removed vs Original", "Reduced vs Removed"))

  # Create comprehensive Jaccard comparison plot
  p_jaccard_comparison <- ggplot(jaccard_long, aes(x = Comparison, y = Jaccard_Index, color = Comparison)) +
    geom_boxplot(position = "dodge", alpha = 0.3, fill = "cornsilk", outlier.shape = NA, size = 0.3) +
    geom_point(aes(shape = as.factor(Iteration)), position = position_jitter(width = 0.1), alpha = 0.7, size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "chartreuse3", size = 0.5) +      # Perfect overlap
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", size = 0.3) +         # Moderate overlap
    geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 0.3, alpha = 0.5) + # No overlap
    facet_grid(Correction_Type ~ ParameterCombination, scales = "free_x") +
    scale_color_colorblind() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
    theme_light() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11)
    ) +
    labs(
      x = "Comparison",
      y = "Jaccard Index",
      title = "Feature Selection Overlap - All Parameter Combinations",
      subtitle = paste(format(nTrials, scientific = FALSE), "trials per iteration,", nSamples, "iterations")
    ) +
    guides(shape = "none", color = guide_legend(title = "Comparison Type"))

  # Combine the two plots side by side for comparison
  p_combined_correlations <- cowplot::plot_grid(
    p_single_trial_correlations,
    p_all_correlations,
    p_jaccard_comparison,
    labels = c("A", "B", "C"),
    rel_widths = c(1, 4, 4),
    align = "h", axis = "b",
    ncol = 3
  )

  # Display and save the combined comparative plot
  print(p_combined_correlations)

  # Save with descriptive filename
  ggsave(paste0("combined_correlations_comparison_",
                format(nTrials, scientific = FALSE), "trials_",
                nSamples, "iterations.svg"),
         p_combined_correlations, width = 12, height = 12)

  # ===============================================================================
  # PARAMETER IMPACT ANALYSIS
  # ===============================================================================

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
      # Extract individual parameters from string like "_OBFALSE_NNSFALSE"
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
}

# Final completion message
cat("\n=== ALL PROCESSING COMPLETE ===\n")
cat("Results stored in workspace variables for further analysis.\n")
cat("Parameter impact analysis results available in 'parameter_impact_analysis' variable.\n")
print(parameter_impact_analysis)