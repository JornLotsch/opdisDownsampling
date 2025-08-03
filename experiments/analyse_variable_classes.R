# ===============================================================================
# KENDALL'S TAU CORRELATION ANALYSIS FOR P-VALUES
# ===============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(cowplot)

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
# DATA EXTRACTION AND PREPARATION
# ===============================================================================

# Extract p-values for 1-trial experiment
pval_orig <- results_experiment_ascending_significance_1Trial$p_values_orig$p_values

pvals_reduced_1Trial <- do.call(cbind.data.frame, lapply(results_experiment_ascending_significance_1Trial$pval_reduced_list, function(x) x$p_values))
pvals_removed_1Trial <- do.call(cbind.data.frame, lapply(results_experiment_ascending_significance_1Trial$pval_removed_list, function(x) x$p_values))


# Function for extracting p-values from n-trials experiment
extract_pval_reduced_list_to_df <- function(pval_reduced_list) {
  if (is.null(pval_reduced_list) || length(pval_reduced_list) == 0) {
    return(NULL)
  }

  # Extract all p_values vectors named by iterations
  pval_vectors <- lapply(pval_reduced_list, function(x) x[["p_values"]])

  # Remove NULLs if any iteration missing p_values
  pval_vectors <- pval_vectors[!sapply(pval_vectors, is.null)]

  if (length(pval_vectors) == 0) {
    return(NULL)
  }

  # Ensure consistent vector lengths by padding with NA
  max_len <- max(sapply(pval_vectors, length))
  padded_vectors <- lapply(pval_vectors, function(vec) {
    length(vec) <- max_len # Pads with NA if vec shorter than max_len
    vec
  })

  df <- as.data.frame(padded_vectors)
  names(df) <- names(pval_vectors)

  return(df)
}

# Extract p-values for n-trials experiment
pvals_reduced_nTrials <- lapply(results_experiment_ascending_significance_nTrials, function(group) {
  pval_reduced_list <- group[["pval_reduced_list"]]
  extract_pval_reduced_list_to_df(pval_reduced_list)
})

pvals_removed_nTrials <- lapply(results_experiment_ascending_significance_nTrials, function(group) {
  pval_reduced_list <- group[["pval_removed_list"]]
  extract_pval_reduced_list_to_df(pval_reduced_list)
})

# Remove NULL entries
pvals_reduced_nTrials <- pvals_reduced_nTrials[!sapply(pvals_reduced_nTrials, is.null)]
pvals_removed_nTrials <- pvals_removed_nTrials[!sapply(pvals_removed_nTrials, is.null)]

# Create variable class labels
var_class <- gsub("[^A-Za-z]", "", names(actual_data))[gsub("[^A-Za-z]", "", names(actual_data)) != gsub("[^A-Za-z]", "", actual_class)]

# Remove C type features 
# Identify significant features for original data using p-values
if (use_ABC_for_feature_selection && requireNamespace("ABCanalysis", quietly = TRUE)) {
  # Use ABCanalysis for feature selection based on -log10 transformed p-values
  feature_indices_orig_C <- safe_test_execution(
    function() ABCanalysis::ABCanalysis(-log10(pmax(pval_orig, 1e-100)))$Aind,
    "ABCanalysis for raw p-values",
    which(pval_orig < significance_threshold)
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
      pvals_reduced_1Trial <- pvals_reduced_1Trial[-feature_indices_orig_C,, drop = FALSE]
      pvals_removed_1Trial <- pvals_removed_1Trial[-feature_indices_orig_C,, drop = FALSE]

      pvals_reduced_nTrials <- lapply(pvals_reduced_nTrials, function(pvals_reduced_nTrials) pvals_reduced_nTrials[-feature_indices_orig_C,, drop = FALSE])
      pvals_removed_nTrials <- lapply(pvals_removed_nTrials, function(pvals_removed_nTrials) pvals_removed_nTrials[-feature_indices_orig_C,, drop = FALSE])

      var_class <- var_class[-feature_indices_orig_C]

      cat(paste("ABC filtering removed", length(feature_indices_orig_C), "features\n"))
    }
  }
}


# ===============================================================================
# KENDALL'S TAU CALCULATION FUNCTIONS
# ===============================================================================

#' Calculate Kendall's tau correlation between original and downsampled p-values by variable class
#' @param pval_orig Original p-values vector
#' @param pvals_downsampled Data frame of downsampled p-values
#' @param var_class Variable class labels
#' @param dataset_name Name for the dataset (e.g., "Reduced", "Removed")
#' @return Data frame with Kendall's tau values by variable class
calculate_kendall_by_class <- function(pval_orig, pvals_downsampled, var_class, dataset_name) {
  if (is.null(pvals_downsampled) || ncol(pvals_downsampled) == 0) {
    return(NULL)
  }

  # Combine data for analysis
  df_combined <- cbind.data.frame(pval_orig = pval_orig, pvals_downsampled, var_class = var_class)

  # Calculate Kendall's tau by variable class
  kendall_results <- df_combined %>%
    group_by(var_class) %>%
    summarise(across(all_of(names(pvals_downsampled)), ~ cor(pval_orig, ., method = "kendall", use = "complete.obs")),
              .groups = "drop")

  kendall_results$DataSet <- dataset_name
  return(as.data.frame(kendall_results))
}

# ===============================================================================
# PROCESS 1-TRIAL EXPERIMENT
# ===============================================================================

# Calculate Kendall's tau for 1-trial experiment
kendall_reduced_1Trial <- calculate_kendall_by_class(pval_orig, pvals_reduced_1Trial, var_class, "Reduced")
kendall_removed_1Trial <- calculate_kendall_by_class(pval_orig, pvals_removed_1Trial, var_class, "Removed")

# Combine 1-trial results
kendall_by_class_1Trial <- rbind.data.frame(kendall_reduced_1Trial, kendall_removed_1Trial)
kendall_by_class_1Trial$Parameters <- "OpF_NoF"
kendall_by_class_1Trial$SingleMultiple <- "Single"

# ===============================================================================
# PROCESS N-TRIALS EXPERIMENT
# ===============================================================================

# Calculate Kendall's tau for n-trials experiment - Reduced
kendall_reduced_nTrials <- lapply(names(pvals_reduced_nTrials), function(param_name) {
  result <- calculate_kendall_by_class(pval_orig, pvals_reduced_nTrials[[param_name]], var_class, "Reduced")
  if (!is.null(result)) {
    result$Parameters <- param_name
  }
  return(result)
})
kendall_reduced_nTrials <- do.call(rbind, kendall_reduced_nTrials[!sapply(kendall_reduced_nTrials, is.null)])

# Calculate Kendall's tau for n-trials experiment - Removed
kendall_removed_nTrials <- lapply(names(pvals_removed_nTrials), function(param_name) {
  result <- calculate_kendall_by_class(pval_orig, pvals_removed_nTrials[[param_name]], var_class, "Removed")
  if (!is.null(result)) {
    result$Parameters <- param_name
  }
  return(result)
})
kendall_removed_nTrials <- do.call(rbind, kendall_removed_nTrials[!sapply(kendall_removed_nTrials, is.null)])

# Combine n-trials results
kendall_by_class_nTrials <- rbind.data.frame(kendall_reduced_nTrials, kendall_removed_nTrials)
kendall_by_class_nTrials$SingleMultiple <- "Multiple"

# ===============================================================================
# COMBINE ALL RESULTS
# ===============================================================================

kendall_by_class_all <- rbind.data.frame(kendall_by_class_1Trial, kendall_by_class_nTrials)

# Convert to long format for plotting
kendall_by_class_all_long <- reshape2::melt(kendall_by_class_all,
                                            id.vars = c("var_class", "DataSet", "Parameters", "SingleMultiple"),
                                            variable.name = "Iteration",
                                            value.name = "KendallTau")

# ===============================================================================
# CREATE IMPROVED PLOTS
# ===============================================================================

# Define consistent colors using colorblind palette
dataset_colors <- scale_color_colorblind()$palette(2)
names(dataset_colors) <- c("Reduced", "Removed")

# Calculate mode function
calculate_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Common theme matching your style
common_theme <- theme_light() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0.5)),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(colour = "black"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9),
    plot.title = element_text(hjust = 0)
  )

# Function to add annotations if requested
add_kendall_annotations <- function(p, data, add_annotations = TRUE, rotate_text = FALSE) {
  if (!add_annotations) return(p)

  # For multiple plot with facets, group by Parameters as well
  if ("Parameters" %in% names(data)) {
    # Calculate statistics for each combination of var_class, DataSet, AND Parameters
    stats_by_comparison <- data %>%
      group_by(var_class, DataSet, Parameters) %>%
      summarise(
        median_val = median(KendallTau, na.rm = TRUE),
        variance_val = var(KendallTau, na.rm = TRUE),
        min_val = min(KendallTau, na.rm = TRUE),
        mode_val = calculate_mode(KendallTau),
        .groups = "drop"
      )
  } else {
    # For single plot, only group by var_class and DataSet
    stats_by_comparison <- data %>%
      group_by(var_class, DataSet) %>%
      summarise(
        median_val = median(KendallTau, na.rm = TRUE),
        variance_val = var(KendallTau, na.rm = TRUE),
        min_val = min(KendallTau, na.rm = TRUE),
        mode_val = calculate_mode(KendallTau),
        .groups = "drop"
      )
  }

  # Create annotation data frame
  annotation_df <- data.frame()

  # Loop through each row of statistics
  for (i in 1:nrow(stats_by_comparison)) {
    stats <- stats_by_comparison[i,]

    # Position annotations above boxes
    var_class_num <- which(unique(data$var_class) == stats$var_class)
    x_offset <- ifelse(stats$DataSet == "Reduced", -0.2, 0.2)

    comp_annotations <- data.frame(
      var_class = rep(stats$var_class, 4),
      DataSet = rep(stats$DataSet, 4),
      x = rep(var_class_num + x_offset, 4),
      y = c(1.2, 1.3, 1.4, 1.5), # Stacked vertically
      label = c(
        paste("Med:", round(stats$median_val, 3)),
        paste("Var:", round(stats$variance_val, 3)),
        paste("Min:", round(stats$min_val, 3)),
        paste("Mode:", round(stats$mode_val, 3))
      )
    )

    # Add Parameters column if it exists in the original data
    if ("Parameters" %in% names(data)) {
      comp_annotations$Parameters <- rep(stats$Parameters, 4)
    }

    annotation_df <- rbind(annotation_df, comp_annotations)
  }

  # Add annotations to plot
  if (nrow(annotation_df) > 0) {
    p <- p +
      geom_text(
        data = annotation_df,
        aes(x = x, y = y, label = label, color = DataSet),
        inherit.aes = FALSE,
        hjust = 0.5, vjust = 0.5, size = 1.5, fontface = "plain",
        angle = ifelse(rotate_text, 90, 0)
      ) +
      scale_y_continuous(
        limits = c(min(kendall_by_class_all_long$KendallTau, na.rm = TRUE) - 0.1, 1.55),
        breaks = seq(from = floor(min(kendall_by_class_all_long$KendallTau, na.rm = TRUE) * 4) / 4,
                     to = 1, by = 0.25)
      )
  }

  return(p)
}


# Create plot for Single trial (NO shape legend, with optional annotations)
p_single <- kendall_by_class_all_long %>%
  filter(SingleMultiple == "Single") %>%
  ggplot(aes(x = var_class, y = KendallTau, color = DataSet)) +
  geom_boxplot(position = position_dodge(width = 0.7), alpha = 0.3,
               outlier.shape = NA, size = 0.4, fill = "white") +
  geom_point(aes(shape = Iteration),
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2),
             alpha = 0.8, size = 1.2) +
  scale_color_colorblind() +
  scale_shape_manual(values = 0:19) +
  scale_y_continuous(limits = c(min(kendall_by_class_all_long$KendallTau, na.rm = TRUE), 1),
                     breaks = seq(from = floor(min(kendall_by_class_all_long$KendallTau, na.rm = TRUE) * 4) / 4,
                                  to = 1, by = 0.25)) +
  common_theme +
  labs(
    x = "Variable Class",
    y = "Kendall's τ",
    title = "Single Trial",
    subtitle = "Kendall's τ between original and downsampled p-values"
  ) +
  guides(shape = "none") # Remove shape legend

# Add annotations to single plot
p_single <- add_kendall_annotations(p_single,
                                    kendall_by_class_all_long %>% filter(SingleMultiple == "Single"),
                                    add_annotations = add_annotations,
                                    rotate_text = TRUE)

# Create plot for Multiple trials (WITH shape legend, with optional annotations)
p_multiple <- kendall_by_class_all_long %>%
  filter(SingleMultiple == "Multiple") %>%
  ggplot(aes(x = var_class, y = KendallTau, color = DataSet)) +
  geom_boxplot(position = position_dodge(width = 0.7), alpha = 0.3,
               outlier.shape = NA, size = 0.4, fill = "white") +
  geom_point(aes(shape = Iteration),
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.15),
             alpha = 0.7, size = 0.9) +
  facet_wrap(~Parameters, ncol = 2, scales = "free_x") +
  scale_color_colorblind() +
  scale_shape_manual(values = 0:19) +
  scale_y_continuous(limits = c(min(kendall_by_class_all_long$KendallTau, na.rm = TRUE), 1.55),
                     breaks = seq(from = floor(min(kendall_by_class_all_long$KendallTau, na.rm = TRUE) * 4) / 4,
                                  to = 1, by = 0.25)) +
  common_theme +
  labs(
    x = "Variable Class",
    y = "Kendall's τ",
    title = "P-Value Rank Correlation: Multiple Trials",
    subtitle = "Kendall's τ across different parameter combinations",
    shape = "Iteration"
  )

# Add annotations to multiple plot
p_multiple <- add_kendall_annotations(p_multiple,
                                      kendall_by_class_all_long %>% filter(SingleMultiple == "Multiple"),
                                      add_annotations = add_annotations,
                                      rotate_text = FALSE)

# Combine plots with 1 row and rel_widths
p_combined <- plot_grid(
  p_single, p_multiple,
  labels = "AUTO",
  ncol = 2,
  rel_widths = c(1, 4),
  align = "h"
)

# Add overall title with left alignment
# title <- ggdraw() +
#   draw_label(
#     "Kendall's Tau Correlation Analysis: P-Value Preservation in Downsampling",
#     fontface = 'bold',
#     size = 14,
#     hjust = 0,
#     x = 0.02
#   )

# Final combined plot
# p_final <- plot_grid(
#   title, p_combined,
#   ncol = 1,
#   rel_heights = c(0.1, 1)
# )
# ===============================================================================
# DISPLAY RESULTS
# ===============================================================================

print(p_combined)

# ===============================================================================
# SUMMARY STATISTICS
# ===============================================================================

# Calculate summary statistics
summary_stats <- kendall_by_class_all_long %>%
  group_by(var_class, DataSet, SingleMultiple, Parameters) %>%
  summarise(
    mean_tau = mean(KendallTau, na.rm = TRUE),
    median_tau = median(KendallTau, na.rm = TRUE),
    sd_tau = sd(KendallTau, na.rm = TRUE),
    min_tau = min(KendallTau, na.rm = TRUE),
    max_tau = max(KendallTau, na.rm = TRUE),
    n_obs = sum(!is.na(KendallTau)),
    .groups = "drop"
  )

cat("\nSummary Statistics for Kendall's Tau:\n")
print(summary_stats)

# Export plots and data if needed
# ggsave("kendall_tau_analysis.png", p_final, width = 12, height = 10, dpi = 300)
# write.csv(kendall_by_class_all_long, "kendall_tau_data.csv", row.names = FALSE)

cat("\nKendall's Tau correlation analysis completed successfully.\n")

# Save the plot
output_filename <- paste0("Tau_by_variable_class_",
                          format(nTrials, scientific = FALSE), "trials_",
                          nSamples, "iterations_", downsampling_size, "sampled",
                          "_ABC_for_feature_selection", use_ABC_for_feature_selection, ".svg")

ggsave(output_filename, p_combined, width = 14, height = 14)
cat(sprintf("Tau_by_variable_class plot saved as %s", output_filename))