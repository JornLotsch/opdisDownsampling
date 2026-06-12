#' Variable Selection Functions
#'
#' Helper function for selecting variables using Principal Component Analysis.
#'
#' This function optionally performs PCA-based variable selection. If
#' \code{PCAimportance = TRUE}, variables are ranked by the absolute loadings in
#' the first principal components and variables above the median importance score
#' are selected. If PCA is disabled or fails, all variables are used.
#'
#' @param DataAndClasses A data frame containing the data and class labels.
#' @param PCAimportance A logical value indicating whether to use PCA to identify
#'   relevant variables.
#' @param list.of.seeds A vector of integer values used for reproducibility.
#'
#' @return A character vector of selected variable names.
#'
#' @details The function excludes the class column from variable selection. If no
#' variables are selected, it falls back to using all variables.
#'
#' @importFrom stats prcomp
select_variables <- function(DataAndClasses, PCAimportance,
                             list.of.seeds) {
  # Get all variable names (excluding class column)
  all_vars <- names(DataAndClasses)[1:(ncol(DataAndClasses) - 1)]

  # PCA-based variable selection
  if (PCAimportance) {
    # Extract numeric data for PCA (excluding class column)
    numeric_data <- DataAndClasses[, all_vars, drop = FALSE]

    # Perform PCA
    pca_result <- tryCatch(
      {
        prcomp(numeric_data, scale. = TRUE, center = TRUE)
      },
      error = function(e) {
        warning("PCA failed, using all variables: ", e$message)
        NULL
      }
    )

    if (!is.null(pca_result)) {
      # Select variables based on PCA loadings
      # Use variables with high absolute loadings in first few components
      n_components <- min(3, ncol(pca_result$rotation)) # Use up to 3 components

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

  # Ensure we have at least some variables selected
  if (length(selectedVars) == 0) {
    warning("No variables selected by criteria, using all variables")
    selectedVars <- all_vars
  }

  return(selectedVars)
}
