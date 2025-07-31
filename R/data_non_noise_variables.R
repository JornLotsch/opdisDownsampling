#' Simple noise detection using KS and Box tests only
#'
#' @param var_data Numeric vector of variable data
#' @param test_stat Statistical test to use for uniformity (currently only "ks" implemented)
#' @param alpha Significance threshold (default: 0.05)
#' @param seed Random seed for reproducible tests (optional)
#' @return TRUE if variable is likely noise, FALSE if signal
#'
#' @importFrom stats ks.test runif median Box.test
#'
simple_noise_detection <- function(var_data, test_stat = "ks", alpha = 0.05, seed = NULL) {
  # Set seed if provided for reproducible uniform reference
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Remove NAs and check data quality
  var_data <- var_data[!is.na(var_data)]

  if (length(var_data) < 21) {
    return(FALSE)  # Insufficient data - treat as signal
  }

  # Get range
  min_val <- min(var_data)
  max_val <- max(var_data)

  # Check for constant variable
  if (min_val == max_val) {
    return(FALSE)  # Constant variable - treat as signal
  }

  tryCatch({
    # Test 1: Uniformity test (KS test against uniform distribution)
    uniform_test <- ks.test(var_data, "punif", min_val, max_val)
    uniformity_p <- uniform_test$p.value

    # Test 2: Independence test (Ljung-Box for autocorrelation)
    if (length(var_data) > 10) {
      lb_test <- Box.test(var_data, lag = min(10, length(var_data) - 1), type = "Ljung-Box")
      independence_p <- lb_test$p.value
    } else {
      independence_p <- 1.0  # Assume independent for very small samples
    }

    # Variable is noise if it looks uniform AND shows no autocorrelation
    # Use a more lenient threshold for uniformity (twice the alpha)
    uniformity_pass <- uniformity_p >= (alpha * 2)  # More lenient
    independence_pass <- independence_p > alpha

    # Both tests must suggest noise
    is_noise <- uniformity_pass && independence_pass

    return(is_noise)

  }, error = function(e) {
    return(FALSE)  # Treat errors as signal variables
  })
}

#' Identify signal variables using simple noise detection
#'
#' @param data_df Data frame with variables to test
#' @param test_stat Statistical test for uniformity test (default: "ks")
#' @param significance_threshold Threshold for noise detection (default: 0.05)
#' @param verbose Print diagnostic information
#' @param seed Random seed for reproducible tests (optional)
#'
#' @return Vector of variable names that are likely signal (not noise)
identify_non_noise_variables <- function(data_df,
                                           test_stat = "ks",
                                           significance_threshold = 0.05,
                                           verbose = FALSE,
                                           seed = NULL) {

  # Extract numeric variables (exclude class column if present)
  if ("Cls" %in% names(data_df)) {
    numeric_vars <- data_df[, !names(data_df) %in% "Cls", drop = FALSE]
  } else {
    numeric_vars <- data_df
  }

  # Apply simple noise detection with seed forwarding
  is_noise <- apply(numeric_vars, 2, function(x) {
    simple_noise_detection(x, test_stat, significance_threshold, seed)
  })

  # Signal variables are those that are NOT noise
  signal_vars <- names(is_noise)[!is_noise]

  if (verbose) {
    cat(paste("Simple noise detection (KS + Ljung-Box):\n"))
    cat("Signal variables (selected):", length(signal_vars), "out of", length(is_noise), "\n")
    cat("Noise variables (filtered):", sum(is_noise), "out of", length(is_noise), "\n")

    if (length(signal_vars) == 0) {
      cat("Warning: No signal variables detected. Using all variables.\n")
    }
  }

  return(signal_vars)
}