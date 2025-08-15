#' Variable Selection Functions
#'
#' Helper functions for selecting relevant variables using PCA and/or
#' non-uniform distribution tests.

#' Perform Non-uniform Variable Selection
#'
#' Internal function that performs variable selection based on distribution uniformity tests.
#'
#' @param DataAndClasses A data frame containing the data and class labels
#' @param selectedVars Currently selected variables (from PCA if applicable)
#' @param PCAimportance Whether PCA was used for initial selection
#' @param NonNoiseSelection Whether to perform non-uniform selection
#' @param UniformTestStat Statistical test for uniformity
#' @param UniformThreshold Threshold for uniformity test
#' @param list.of.seeds Vector of seeds for reproducibility
#'
#' @return Vector of selected variable names
perform_nonuniform_selection <- function(DataAndClasses, selectedVars, PCAimportance,
                                         NonNoiseSelection, UniformTestStat, UniformThreshold,
                                         list.of.seeds) {
  if (!NonNoiseSelection) {
    return(selectedVars)
  }

  # Get all available variables (excluding class column)
  all_vars <- names(DataAndClasses)[1:(ncol(DataAndClasses) - 1)]

  # Start with selectedVars or all variables if no PCA selection
  vars_to_test <- if (PCAimportance) selectedVars else all_vars

  # Perform uniformity test for each variable
  uniform_p_values <- sapply(vars_to_test, function(var_name) {
    var_data <- DataAndClasses[[var_name]]

    # Test for uniformity - lower p-value means less uniform (more informative)
    tryCatch({
      if (UniformTestStat == "ks") {
        # Kolmogorov-Smirnov test against uniform distribution
        ks.test(var_data, "punif", min = min(var_data), max = max(var_data))$p.value
      } else if (UniformTestStat == "ad") {
        # Anderson-Darling test - would need additional implementation
        # For now, fallback to KS test
        ks.test(var_data, "punif", min = min(var_data), max = max(var_data))$p.value
      } else {
        # Default to KS test for unsupported test statistics
        ks.test(var_data, "punif", min = min(var_data), max = max(var_data))$p.value
      }
    }, error = function(e) {
      # Return high p-value (uniform) if test fails
      1.0
    })
  })

  # Select variables with p-value below threshold (non-uniform = informative)
  selected_by_uniformity <- names(uniform_p_values)[uniform_p_values < UniformThreshold]

  return(selected_by_uniformity)
}

#' Select Variables Using PCA and/or Non-uniform Distribution Tests
#'
#' This function performs variable selection using Principal Component Analysis (PCA)
#' and/or non-uniform distribution tests to identify relevant variables for analysis.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param PCAimportance A logical value indicating whether to use PCA to identify
#'   relevant variables.
#' @param NonNoiseSelection A logical value indicating whether to use non-uniform
#'   distribution tests to identify relevant variables.
#' @param UniformTestStat A character string specifying the statistical test to be used for
#'   non-uniform variable selection.
#' @param UniformThreshold Threshold value for non-uniform variable selection.
#' @param list.of.seeds A vector of integer values used for reproducibility.
#'
#' @return A character vector of selected variable names.
#'
#' @details The function operates in two phases:
#' \itemize{
#'   \item PCA Selection: If enabled, performs PCA and selects variables based on loadings
#'   \item Non-uniform Selection: If enabled, tests variables for uniformity and selects non-uniform ones
#' }
#'
#' @importFrom stats prcomp
select_variables <- function(DataAndClasses, PCAimportance, NonNoiseSelection,
                             UniformTestStat, UniformThreshold, list.of.seeds) {

  # Get all variable names (excluding class column)
  all_vars <- names(DataAndClasses)[1:(ncol(DataAndClasses) - 1)]

  # Phase 1: PCA-based variable selection
  if (PCAimportance) {
    # Extract numeric data for PCA (excluding class column)
    numeric_data <- DataAndClasses[, all_vars, drop = FALSE]

    # Perform PCA
    pca_result <- tryCatch({
      prcomp(numeric_data, scale. = TRUE, center = TRUE)
    }, error = function(e) {
      warning("PCA failed, using all variables: ", e$message)
      NULL
    })

    if (!is.null(pca_result)) {
      # Select variables based on PCA loadings
      # Use variables with high absolute loadings in first few components
      n_components <- min(3, ncol(pca_result$rotation))  # Use up to 3 components

      # Calculate importance score for each variable
      var_importance <- rowSums(abs(pca_result$rotation[, 1:n_components, drop = FALSE]))

      # Select top variables (e.g., above median importance)
      importance_threshold <- median(var_importance)
      selectedVars <- names(var_importance)[var_importance > importance_threshold]

    } else {
      # Fallback if PCA fails
      selectedVars <- all_vars
    }
  } else {
    # Start with all variables if no PCA selection
    selectedVars <- all_vars
  }

  # Phase 2: Non-uniform distribution selection
  selectedVars <- perform_nonuniform_selection(
    DataAndClasses = DataAndClasses,
    selectedVars = selectedVars,
    PCAimportance = PCAimportance,
    NonNoiseSelection = NonNoiseSelection,
    UniformTestStat = UniformTestStat,
    UniformThreshold = UniformThreshold,
    list.of.seeds = list.of.seeds
  )

  # Ensure we have at least some variables selected
  if (length(selectedVars) == 0) {
    warning("No variables selected by criteria, using all variables")
    selectedVars <- all_vars
  }

  return(selectedVars)
}