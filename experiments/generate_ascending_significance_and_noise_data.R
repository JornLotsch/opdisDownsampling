# ---------------------------------------------------------------------------
# DATA GENERATION FUNCTIONS AND PARAMETERS
# ---------------------------------------------------------------------------

# Default parameters for synthetic data generation
DEFAULT_RANDOM_SEED <- 42 # Random seed for reproducibility
DEFAULT_SAMPLES_PER_CLASS <- 100 # Number of samples per class
DEFAULT_NUM_VARIABLES <- 50 # Number of variables (features) per data type
DEFAULT_SD <- 5 # Standard deviation for Gaussian variables
DEFAULT_BASE_RANGE <- 100 # Base range for uniform distribution
DEFAULT_CLASS_WIDTH_RATIO <- 0.2 # Width ratio for uniform class distributions

#' Generate synthetic dataset with ascending class separation
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_range Numeric vector of length 2 defining the range for class 1 means
#' @param class2_range Numeric vector of length 2 defining the range for class 2 means
#' @param sd Standard deviation for normal distributions (default: 5)
#'
#' @return A data frame with class labels and generated variables
generate_synthetic_data <- function(n_samples, n_vars, seed,
                                    class1_range = c(1, 50),
                                    class2_range = c(1, 48),
                                    sd = DEFAULT_SD) {

  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (length(class1_range) != 2 || length(class2_range) != 2) {
    stop("class1_range and class2_range must be vectors of length 2")
  }

  # Pre-calculate means to avoid repeated computation
  class1_means <- seq(class1_range[1], class1_range[2], length.out = n_vars)
  class2_means <- seq(class2_range[1], class2_range[2], length.out = n_vars)

  # Generate variables with unique seeds for each
  variables <- lapply(1:n_vars, function(i) {
    # Use different seed for each variable to ensure independence
    set.seed(seed)
    class1_data <- rnorm(n_samples, mean = class1_means[i], sd = sd)
    class2_data <- rnorm(n_samples, mean = class2_means[i], sd = sd)
    c(class1_data, class2_data)
  })

  # Create data frame with proper column names
  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("Var", 1:n_vars)

  data.frame(
    Cls = rep(1:2, each = n_samples),
    var_matrix
  )
}

#' Creates a dataset with Gaussian mixture, noise, and uniform variables
#'
#' @param random_seed Random seed for reproducibility (default: DEFAULT_RANDOM_SEED)
#' @param samples_per_class Number of samples per class (default: DEFAULT_SAMPLES_PER_CLASS)
#' @param num_variables Number of variables for each type (default: DEFAULT_NUM_VARIABLES)
#' @param sd Standard deviation for Gaussian variables (default: DEFAULT_SD)
#' @param base_range Base range for uniform distributions (default: DEFAULT_BASE_RANGE)
#' @param class_width_ratio Width ratio for uniform class distributions (default: DEFAULT_CLASS_WIDTH_RATIO)
#' @param min_class_diff_ratio Minimum class difference as ratio of base_range (default: 0.0005)
#' @param max_class_diff_ratio Maximum class difference as ratio of base_range (default: 0.025)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return A data frame with class labels and three types of variables
generate_comprehensive_synthetic_data <- function(random_seed = DEFAULT_RANDOM_SEED,
                                                  samples_per_class = DEFAULT_SAMPLES_PER_CLASS,
                                                  num_variables = DEFAULT_NUM_VARIABLES,
                                                  sd = DEFAULT_SD,
                                                  base_range = DEFAULT_BASE_RANGE,
                                                  class_width_ratio = DEFAULT_CLASS_WIDTH_RATIO,
                                                  min_class_diff_ratio = 0.0005,
                                                  max_class_diff_ratio = 0.025,
                                                  verbose = TRUE) {

  # Input validation
  if (samples_per_class <= 0 || num_variables <= 0) {
    stop("samples_per_class and num_variables must be positive integers")
  }
  if (sd <= 0 || base_range <= 0) {
    stop("sd and base_range must be positive values")
  }
  if (class_width_ratio <= 0 || class_width_ratio >= 1) {
    stop("class_width_ratio must be between 0 and 1")
  }

  if (verbose) {
    cat("Generating comprehensive synthetic dataset...\n")
    cat(sprintf("- %d samples per class, %d variables per type\n", samples_per_class, num_variables))
  }

  # Generate base Gaussian mixture data with ascending class separation
  if (verbose) cat("- Generating Gaussian mixture variables...\n")
  base_data <- generate_synthetic_data(samples_per_class, num_variables, random_seed)

  # Pre-calculate ranges for noise variables to avoid repeated computation
  class1_means <- seq(1, 50, length.out = num_variables)
  class2_means <- seq(1, 48, length.out = num_variables)
  overall_means <- (class1_means + class2_means) / 2
  range_mins <- overall_means - 3 * sd
  range_maxs <- overall_means + 3 * sd

  # Generate pure noise variables (uniform random within same value range as informative features)
  if (verbose) cat("- Generating noise variables...\n")
  noise_data <- lapply(1:num_variables, function(i) {
    # Use systematic seed generation to avoid collisions
    set.seed(random_seed + 10000 + i)
    runif(2 * samples_per_class, min = range_mins[i], max = range_maxs[i])
  })

  # Pre-calculate class differences for uniform variables
  class_differences <- seq(min_class_diff_ratio, max_class_diff_ratio,
                           length.out = num_variables) * base_range
  class1_center <- 50 # Fixed center for class 1
  class2_centers <- class1_center + class_differences
  class_width <- base_range * class_width_ratio

  # Generate uniform variables with very subtle but ascending significant differences between classes
  if (verbose) cat("- Generating uniform variables with class differences...\n")
  uniform_class_diff_data <- lapply(1:num_variables, function(i) {
    # Use systematic seed generation to avoid collisions
    set.seed(random_seed + 20000 + i)

    class1_data <- runif(samples_per_class,
                         min = class1_center - class_width / 2,
                         max = class1_center + class_width / 2)

    class2_data <- runif(samples_per_class,
                         min = class2_centers[i] - class_width / 2,
                         max = class2_centers[i] + class_width / 2)

    c(class1_data, class2_data)
  })

  # Combine all three data types into final dataset
  if (verbose) cat("- Combining datasets...\n")

  # Create matrices with proper column names
  noise_matrix <- do.call(cbind, noise_data)
  colnames(noise_matrix) <- paste0("NoiseVar", 1:num_variables)

  uniform_matrix <- do.call(cbind, uniform_class_diff_data)
  colnames(uniform_matrix) <- paste0("UniformVar", 1:num_variables)

  # Combine all data
  final_data <- cbind(base_data, noise_matrix, uniform_matrix)

  # Update Gaussian variable names for clarity
  colnames(final_data)[2:(num_variables + 1)] <- paste0("GMMVar", 1:num_variables)

  if (verbose) {
    cat(sprintf("- Dataset generated: %d observations, %d variables total\n",
                nrow(final_data), ncol(final_data) - 1))
    cat("- Variable types: GMMVar (Gaussian), NoiseVar (Pure noise), UniformVar (Class differences)\n")
  }

  return(final_data)
}

#' Generate dataset with custom variable type proportions
#'
#' @param random_seed Random seed for reproducibility
#' @param samples_per_class Number of samples per class
#' @param n_gmm Number of Gaussian mixture variables
#' @param n_noise Number of pure noise variables
#' @param n_uniform Number of uniform variables with class differences
#' @param ... Additional parameters passed to generate_comprehensive_synthetic_data
#'
#' @return A data frame with specified proportions of variable types
generate_custom_synthetic_data <- function(random_seed = DEFAULT_RANDOM_SEED,
                                           samples_per_class = DEFAULT_SAMPLES_PER_CLASS,
                                           n_gmm = DEFAULT_NUM_VARIABLES,
                                           n_noise = DEFAULT_NUM_VARIABLES,
                                           n_uniform = DEFAULT_NUM_VARIABLES,
                                           ...) {

  # Generate each type separately with different seeds to ensure independence
  gmm_data <- generate_synthetic_data(samples_per_class, n_gmm, random_seed, ...)

  # Generate noise data
  set.seed(random_seed + 100000)
  noise_vars <- lapply(1:n_noise, function(i) {
    set.seed(random_seed + 100000 + i)
    runif(2 * samples_per_class, min = -30, max = 80) # Wide range covering GMM range
  })
  noise_matrix <- do.call(cbind, noise_vars)
  colnames(noise_matrix) <- paste0("NoiseVar", 1:n_noise)

  # Generate uniform variables with class differences
  uniform_vars <- lapply(1:n_uniform, function(i) {
    set.seed(random_seed + 200000 + i)
    class_diff <- (i / n_uniform) * 5 # Increasing class difference

    class1_data <- runif(samples_per_class, min = 45, max = 55)
    class2_data <- runif(samples_per_class, min = 45 + class_diff, max = 55 + class_diff)

    c(class1_data, class2_data)
  })
  uniform_matrix <- do.call(cbind, uniform_vars)
  colnames(uniform_matrix) <- paste0("UniformVar", 1:n_uniform)

  # Update GMM variable names
  colnames(gmm_data)[2:(n_gmm + 1)] <- paste0("GMMVar", 1:n_gmm)

  # Combine all data
  final_data <- cbind(gmm_data, noise_matrix, uniform_matrix)

  return(final_data)
}