# ---------------------------------------------------------------------------
# DATA GENERATION FUNCTIONS AND PARAMETERS (EXTENDED VERSION)
# ---------------------------------------------------------------------------

# Default parameters for synthetic data generation
DEFAULT_RANDOM_SEED <- 42 # Random seed for reproducibility
DEFAULT_SAMPLES_PER_CLASS <- 100 # Number of samples per class
DEFAULT_NUM_VARIABLES <- 50 # Number of variables (features) per data type
DEFAULT_SD <- 5 # Standard deviation for Gaussian variables
DEFAULT_BASE_RANGE <- 100 # Base range for uniform distribution
DEFAULT_CLASS_WIDTH_RATIO <- 0.2 # Width ratio for uniform class distributions

# ---------------------------------------------------------------------------
# ORIGINAL DISTRIBUTION FUNCTIONS (EXISTING)
# ---------------------------------------------------------------------------

#' Generate synthetic dataset with ascending class separation (Component 1)
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_range Numeric vector of length 2 defining the range for class 1 means
#' @param class2_range Numeric vector of length 2 defining the range for class 2 means
#' @param sd Standard deviation for normal distributions (default: 5)
#'
#' @return A data frame with class labels and generated variables
generate_ascending_separation_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                               n_vars = DEFAULT_NUM_VARIABLES,
                                               seed = DEFAULT_RANDOM_SEED,
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

  # Generate variables with same seed for stepwise increasing statistical difference
  variables <- lapply(1:n_vars, function(i) {
    # Use same seed for all variables to ensure stepwise increasing difference pattern
    set.seed(seed)
    class1_data <- rnorm(n_samples, mean = class1_means[i], sd = sd)
    class2_data <- rnorm(n_samples, mean = class2_means[i], sd = sd)
    c(class1_data, class2_data)
  })

  # Create data frame with proper column names
  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("GMMVar", 1:n_vars)

  data.frame(
    Cls = rep(1:2, each = n_samples),
    var_matrix
  )
}

#' Generate Gaussian noise data without class differences (Component 2)
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param sd Standard deviation for normal distributions (default: 5)
#' @param mean_range Range for random means (default: c(1, 50))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_gaussian_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                         n_vars = DEFAULT_NUM_VARIABLES,
                                         seed = DEFAULT_RANDOM_SEED,
                                         sd = DEFAULT_SD,
                                         mean_range = c(1, 50)) {

  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (sd <= 0) {
    stop("sd must be positive")
  }

  # Generate variables with different seeds to ensure independence
  variables <- lapply(1:n_vars, function(i) {
    # Use different seed for each variable to ensure independence
    set.seed(seed + 10000 + i)

    # Random mean for this variable (same for both classes)
    var_mean <- runif(1, min = mean_range[1], max = mean_range[2])

    # Generate data for both classes with same mean (no class difference)
    class1_data <- rnorm(n_samples, mean = var_mean, sd = sd)
    class2_data <- rnorm(n_samples, mean = var_mean, sd = sd)
    c(class1_data, class2_data)
  })

  # Create matrix with proper column names
  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseGaussVar", 1:n_vars)

  return(var_matrix)
}

#' Generate uniform data with ascending class differences (Component 4)
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param base_range Base range for uniform distributions (default: 100)
#' @param class_width_ratio Width ratio for uniform class distributions (default: 0.2)
#' @param min_class_diff_ratio Minimum class difference as ratio of base_range (default: 0.0005)
#' @param max_class_diff_ratio Maximum class difference as ratio of base_range (default: 0.025)
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_uniform_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                             n_vars = DEFAULT_NUM_VARIABLES,
                                             seed = DEFAULT_RANDOM_SEED,
                                             base_range = DEFAULT_BASE_RANGE,
                                             class_width_ratio = DEFAULT_CLASS_WIDTH_RATIO,
                                             min_class_diff_ratio = 0.0005,
                                             max_class_diff_ratio = 0.025) {

  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (base_range <= 0) {
    stop("base_range must be positive")
  }
  if (class_width_ratio <= 0 || class_width_ratio >= 1) {
    stop("class_width_ratio must be between 0 and 1")
  }

  # Pre-calculate class differences for uniform variables (ascending pattern)
  class_differences <- seq(min_class_diff_ratio, max_class_diff_ratio,
                           length.out = n_vars) * base_range
  class1_center <- 50 # Fixed center for class 1
  class2_centers <- class1_center + class_differences
  class_width <- base_range * class_width_ratio

  # Generate variables with same seed for stepwise increasing difference pattern
  variables <- lapply(1:n_vars, function(i) {
    # Use same seed for all variables to ensure stepwise increasing difference pattern
    set.seed(seed)

    class1_data <- runif(n_samples,
                         min = class1_center - class_width / 2,
                         max = class1_center + class_width / 2)

    class2_data <- runif(n_samples,
                         min = class2_centers[i] - class_width / 2,
                         max = class2_centers[i] + class_width / 2)

    c(class1_data, class2_data)
  })

  # Create matrix with proper column names
  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("UniformVar", 1:n_vars)

  return(var_matrix)
}

#' Generate uniform random noise data regardless of class (Component 3)
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param value_range Range for uniform distribution (default: c(-30, 80))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_uniform_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                        n_vars = DEFAULT_NUM_VARIABLES,
                                        seed = DEFAULT_RANDOM_SEED,
                                        value_range = c(-30, 80)) {

  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (length(value_range) != 2) {
    stop("value_range must be a vector of length 2")
  }

  # Generate variables with different seeds to ensure independence
  variables <- lapply(1:n_vars, function(i) {
    # Use different seed for each variable to ensure independence
    set.seed(seed + 20000 + i)

    # Generate uniform random data for both classes (no class difference)
    runif(2 * n_samples, min = value_range[1], max = value_range[2])
  })

  # Create matrix with proper column names
  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseUnifVar", 1:n_vars)

  return(var_matrix)
}

# ---------------------------------------------------------------------------
# NEW BIOMEDICAL DISTRIBUTION FUNCTIONS
# ---------------------------------------------------------------------------

#' Generate log-normal data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_meanlog_range Range for log-mean of class 1 (default: c(1, 3))
#' @param class2_meanlog_range Range for log-mean of class 2 (default: c(1, 2.8))
#' @param sdlog Standard deviation for log-normal distribution (default: 0.5)
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_lognormal_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                               n_vars = DEFAULT_NUM_VARIABLES,
                                               seed = DEFAULT_RANDOM_SEED,
                                               class1_meanlog_range = c(1, 3),
                                               class2_meanlog_range = c(1, 2.8),
                                               sdlog = 0.5) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (sdlog <= 0) {
    stop("sdlog must be positive")
  }

  # Pre-calculate log-means
  class1_meanlogs <- seq(class1_meanlog_range[1], class1_meanlog_range[2], length.out = n_vars)
  class2_meanlogs <- seq(class2_meanlog_range[1], class2_meanlog_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rlnorm(n_samples, meanlog = class1_meanlogs[i], sdlog = sdlog)
    class2_data <- rlnorm(n_samples, meanlog = class2_meanlogs[i], sdlog = sdlog)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("LogNormVar", 1:n_vars)
  return(var_matrix)
}

#' Generate log-normal noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param meanlog_range Range for random log-means (default: c(1, 3))
#' @param sdlog Standard deviation for log-normal distribution (default: 0.5)
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_lognormal_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                          n_vars = DEFAULT_NUM_VARIABLES,
                                          seed = DEFAULT_RANDOM_SEED,
                                          meanlog_range = c(1, 3),
                                          sdlog = 0.5) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 30000 + i)
    var_meanlog <- runif(1, min = meanlog_range[1], max = meanlog_range[2])
    # Generate same distribution for both classes
    class1_data <- rlnorm(n_samples, meanlog = var_meanlog, sdlog = sdlog)
    class2_data <- rlnorm(n_samples, meanlog = var_meanlog, sdlog = sdlog)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseLogNormVar", 1:n_vars)
  return(var_matrix)
}

#' Generate binomial data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param n_trials Number of trials for binomial distribution (default: 10)
#' @param class1_prob_range Range for success probability of class 1 (default: c(0.3, 0.7))
#' @param class2_prob_range Range for success probability of class 2 (default: c(0.3, 0.65))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_binomial_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                              n_vars = DEFAULT_NUM_VARIABLES,
                                              seed = DEFAULT_RANDOM_SEED,
                                              n_trials = 10,
                                              class1_prob_range = c(0.3, 0.7),
                                              class2_prob_range = c(0.3, 0.65)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (n_trials <= 0) {
    stop("n_trials must be positive")
  }

  # Pre-calculate probabilities
  class1_probs <- seq(class1_prob_range[1], class1_prob_range[2], length.out = n_vars)
  class2_probs <- seq(class2_prob_range[1], class2_prob_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rbinom(n_samples, size = n_trials, prob = class1_probs[i])
    class2_data <- rbinom(n_samples, size = n_trials, prob = class2_probs[i])
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("BinomVar", 1:n_vars)
  return(var_matrix)
}

#' Generate binomial noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param n_trials Number of trials for binomial distribution (default: 10)
#' @param prob_range Range for random success probabilities (default: c(0.2, 0.8))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_binomial_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                         n_vars = DEFAULT_NUM_VARIABLES,
                                         seed = DEFAULT_RANDOM_SEED,
                                         n_trials = 10,
                                         prob_range = c(0.2, 0.8)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 40000 + i)
    var_prob <- runif(1, min = prob_range[1], max = prob_range[2])
    # Generate same distribution for both classes
    class1_data <- rbinom(n_samples, size = n_trials, prob = var_prob)
    class2_data <- rbinom(n_samples, size = n_trials, prob = var_prob)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseBinomVar", 1:n_vars)
  return(var_matrix)
}

#' Generate Poisson data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_lambda_range Range for lambda parameter of class 1 (default: c(2, 20))
#' @param class2_lambda_range Range for lambda parameter of class 2 (default: c(2, 18))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_poisson_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                             n_vars = DEFAULT_NUM_VARIABLES,
                                             seed = DEFAULT_RANDOM_SEED,
                                             class1_lambda_range = c(2, 20),
                                             class2_lambda_range = c(2, 18)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  # Pre-calculate lambda parameters
  class1_lambdas <- seq(class1_lambda_range[1], class1_lambda_range[2], length.out = n_vars)
  class2_lambdas <- seq(class2_lambda_range[1], class2_lambda_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rpois(n_samples, lambda = class1_lambdas[i])
    class2_data <- rpois(n_samples, lambda = class2_lambdas[i])
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("PoissonVar", 1:n_vars)
  return(var_matrix)
}

#' Generate Poisson noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param lambda_range Range for random lambda parameters (default: c(1, 25))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_poisson_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                        n_vars = DEFAULT_NUM_VARIABLES,
                                        seed = DEFAULT_RANDOM_SEED,
                                        lambda_range = c(1, 25)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 50000 + i)
    var_lambda <- runif(1, min = lambda_range[1], max = lambda_range[2])
    # Generate same distribution for both classes
    class1_data <- rpois(n_samples, lambda = var_lambda)
    class2_data <- rpois(n_samples, lambda = var_lambda)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoisePoissonVar", 1:n_vars)
  return(var_matrix)
}

#' Generate exponential data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_rate_range Range for rate parameter of class 1 (default: c(0.5, 2))
#' @param class2_rate_range Range for rate parameter of class 2 (default: c(0.5, 1.8))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_exponential_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                                 n_vars = DEFAULT_NUM_VARIABLES,
                                                 seed = DEFAULT_RANDOM_SEED,
                                                 class1_rate_range = c(0.5, 2),
                                                 class2_rate_range = c(0.5, 1.8)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  # Pre-calculate rate parameters
  class1_rates <- seq(class1_rate_range[1], class1_rate_range[2], length.out = n_vars)
  class2_rates <- seq(class2_rate_range[1], class2_rate_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rexp(n_samples, rate = class1_rates[i])
    class2_data <- rexp(n_samples, rate = class2_rates[i])
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("ExpVar", 1:n_vars)
  return(var_matrix)
}

#' Generate exponential noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param rate_range Range for random rate parameters (default: c(0.2, 3))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_exponential_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                            n_vars = DEFAULT_NUM_VARIABLES,
                                            seed = DEFAULT_RANDOM_SEED,
                                            rate_range = c(0.2, 3)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 60000 + i)
    var_rate <- runif(1, min = rate_range[1], max = rate_range[2])
    # Generate same distribution for both classes
    class1_data <- rexp(n_samples, rate = var_rate)
    class2_data <- rexp(n_samples, rate = var_rate)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseExpVar", 1:n_vars)
  return(var_matrix)
}

#' Generate gamma data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_shape_range Range for shape parameter of class 1 (default: c(1, 5))
#' @param class2_shape_range Range for shape parameter of class 2 (default: c(1, 4.5))
#' @param rate Fixed rate parameter for both classes (default: 1)
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_gamma_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                           n_vars = DEFAULT_NUM_VARIABLES,
                                           seed = DEFAULT_RANDOM_SEED,
                                           class1_shape_range = c(1, 5),
                                           class2_shape_range = c(1, 4.5),
                                           rate = 1) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (rate <= 0) {
    stop("rate must be positive")
  }

  # Pre-calculate shape parameters
  class1_shapes <- seq(class1_shape_range[1], class1_shape_range[2], length.out = n_vars)
  class2_shapes <- seq(class2_shape_range[1], class2_shape_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rgamma(n_samples, shape = class1_shapes[i], rate = rate)
    class2_data <- rgamma(n_samples, shape = class2_shapes[i], rate = rate)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("GammaVar", 1:n_vars)
  return(var_matrix)
}

#' Generate gamma noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param shape_range Range for random shape parameters (default: c(0.5, 6))
#' @param rate Fixed rate parameter (default: 1)
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_gamma_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                      n_vars = DEFAULT_NUM_VARIABLES,
                                      seed = DEFAULT_RANDOM_SEED,
                                      shape_range = c(0.5, 6),
                                      rate = 1) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 70000 + i)
    var_shape <- runif(1, min = shape_range[1], max = shape_range[2])
    # Generate same distribution for both classes
    class1_data <- rgamma(n_samples, shape = var_shape, rate = rate)
    class2_data <- rgamma(n_samples, shape = var_shape, rate = rate)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseGammaVar", 1:n_vars)
  return(var_matrix)
}

#' Generate beta data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_shape1_range Range for shape1 parameter of class 1 (default: c(1, 5))
#' @param class2_shape1_range Range for shape1 parameter of class 2 (default: c(1, 4.5))
#' @param shape2 Fixed shape2 parameter for both classes (default: 2)
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_beta_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                          n_vars = DEFAULT_NUM_VARIABLES,
                                          seed = DEFAULT_RANDOM_SEED,
                                          class1_shape1_range = c(1, 5),
                                          class2_shape1_range = c(1, 4.5),
                                          shape2 = 2) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (shape2 <= 0) {
    stop("shape2 must be positive")
  }

  # Pre-calculate shape1 parameters
  class1_shape1s <- seq(class1_shape1_range[1], class1_shape1_range[2], length.out = n_vars)
  class2_shape1s <- seq(class2_shape1_range[1], class2_shape1_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rbeta(n_samples, shape1 = class1_shape1s[i], shape2 = shape2)
    class2_data <- rbeta(n_samples, shape1 = class2_shape1s[i], shape2 = shape2)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("BetaVar", 1:n_vars)
  return(var_matrix)
}

#' Generate beta noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param shape1_range Range for random shape1 parameters (default: c(0.5, 6))
#' @param shape2_range Range for random shape2 parameters (default: c(0.5, 6))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_beta_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                     n_vars = DEFAULT_NUM_VARIABLES,
                                     seed = DEFAULT_RANDOM_SEED,
                                     shape1_range = c(0.5, 6),
                                     shape2_range = c(0.5, 6)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 80000 + i)
    var_shape1 <- runif(1, min = shape1_range[1], max = shape1_range[2])
    var_shape2 <- runif(1, min = shape2_range[1], max = shape2_range[2])
    # Generate same distribution for both classes
    class1_data <- rbeta(n_samples, shape1 = var_shape1, shape2 = var_shape2)
    class2_data <- rbeta(n_samples, shape1 = var_shape1, shape2 = var_shape2)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseBetaVar", 1:n_vars)
  return(var_matrix)
}

#' Generate chi-squared data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_df_range Range for degrees of freedom of class 1 (default: c(1, 10))
#' @param class2_df_range Range for degrees of freedom of class 2 (default: c(1, 9))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_chisq_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                           n_vars = DEFAULT_NUM_VARIABLES,
                                           seed = DEFAULT_RANDOM_SEED,
                                           class1_df_range = c(1, 10),
                                           class2_df_range = c(1, 9)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  # Pre-calculate degrees of freedom
  class1_dfs <- seq(class1_df_range[1], class1_df_range[2], length.out = n_vars)
  class2_dfs <- seq(class2_df_range[1], class2_df_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rchisq(n_samples, df = class1_dfs[i])
    class2_data <- rchisq(n_samples, df = class2_dfs[i])
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("ChiSqVar", 1:n_vars)
  return(var_matrix)
}

#' Generate chi-squared noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param df_range Range for random degrees of freedom (default: c(1, 12))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_chisq_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                      n_vars = DEFAULT_NUM_VARIABLES,
                                      seed = DEFAULT_RANDOM_SEED,
                                      df_range = c(1, 12)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 90000 + i)
    var_df <- sample(df_range[1]:df_range[2], 1) # Use integer df
    # Generate same distribution for both classes
    class1_data <- rchisq(n_samples, df = var_df)
    class2_data <- rchisq(n_samples, df = var_df)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseChiSqVar", 1:n_vars)
  return(var_matrix)
}

#' Generate Weibull data with ascending class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param class1_shape_range Range for shape parameter of class 1 (default: c(1, 3))
#' @param class2_shape_range Range for shape parameter of class 2 (default: c(1, 2.8))
#' @param scale Fixed scale parameter for both classes (default: 1)
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_weibull_class_diff_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                             n_vars = DEFAULT_NUM_VARIABLES,
                                             seed = DEFAULT_RANDOM_SEED,
                                             class1_shape_range = c(1, 3),
                                             class2_shape_range = c(1, 2.8),
                                             scale = 1) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }
  if (scale <= 0) {
    stop("scale must be positive")
  }

  # Pre-calculate shape parameters
  class1_shapes <- seq(class1_shape_range[1], class1_shape_range[2], length.out = n_vars)
  class2_shapes <- seq(class2_shape_range[1], class2_shape_range[2], length.out = n_vars)

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    class1_data <- rweibull(n_samples, shape = class1_shapes[i], scale = scale)
    class2_data <- rweibull(n_samples, shape = class2_shapes[i], scale = scale)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("WeibullVar", 1:n_vars)
  return(var_matrix)
}

#' Generate Weibull noise data without class differences
#'
#' @param n_samples Number of samples per class
#' @param n_vars Number of variables to generate
#' @param seed Random seed for reproducibility
#' @param shape_range Range for random shape parameters (default: c(0.5, 4))
#' @param scale_range Range for random scale parameters (default: c(0.5, 2))
#'
#' @return A matrix with n_vars columns and 2*n_samples rows
generate_weibull_noise_data <- function(n_samples = DEFAULT_SAMPLES_PER_CLASS,
                                        n_vars = DEFAULT_NUM_VARIABLES,
                                        seed = DEFAULT_RANDOM_SEED,
                                        shape_range = c(0.5, 4),
                                        scale_range = c(0.5, 2)) {
  # Input validation
  if (n_samples <= 0 || n_vars <= 0) {
    stop("n_samples and n_vars must be positive integers")
  }

  variables <- lapply(1:n_vars, function(i) {
    set.seed(seed + 100000 + i)
    var_shape <- runif(1, min = shape_range[1], max = shape_range[2])
    var_scale <- runif(1, min = scale_range[1], max = scale_range[2])
    # Generate same distribution for both classes
    class1_data <- rweibull(n_samples, shape = var_shape, scale = var_scale)
    class2_data <- rweibull(n_samples, shape = var_shape, scale = var_scale)
    c(class1_data, class2_data)
  })

  var_matrix <- do.call(cbind, variables)
  colnames(var_matrix) <- paste0("NoiseWeibullVar", 1:n_vars)
  return(var_matrix)
}

# ---------------------------------------------------------------------------
# COMPREHENSIVE GENERATION FUNCTIONS (UPDATED)
# ---------------------------------------------------------------------------

#' Generate comprehensive synthetic dataset with all original 4 components
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
#' @return A data frame with class labels and four types of variables
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
    cat("Generating comprehensive synthetic dataset with 4 components...\n")
    cat(sprintf("- %d samples per class, %d variables per component\n", samples_per_class, num_variables))
  }

  # Component 1: Gaussian mixture data with ascending class separation
  if (verbose) cat("- Component 1: Generating Gaussian mixture variables with ascending class separation...\n")
  component1_data <- generate_ascending_separation_data(
    n_samples = samples_per_class,
    n_vars = num_variables,
    seed = random_seed,
    sd = sd
  )

  # Component 2: Gaussian noise without class differences
  if (verbose) cat("- Component 2: Generating Gaussian noise variables without class differences...\n")
  component2_matrix <- generate_gaussian_noise_data(
    n_samples = samples_per_class,
    n_vars = num_variables,
    seed = random_seed,
    sd = sd
  )

  # Component 3: Uniform variables with ascending class differences
  if (verbose) cat("- Component 3: Generating uniform variables with ascending class differences...\n")
  component3_matrix <- generate_uniform_class_diff_data(
    n_samples = samples_per_class,
    n_vars = num_variables,
    seed = random_seed,
    base_range = base_range,
    class_width_ratio = class_width_ratio,
    min_class_diff_ratio = min_class_diff_ratio,
    max_class_diff_ratio = max_class_diff_ratio
  )

  # Component 4: Uniform random noise regardless of class
  if (verbose) cat("- Component 4: Generating uniform random noise variables...\n")
  component4_matrix <- generate_uniform_noise_data(
    n_samples = samples_per_class,
    n_vars = num_variables,
    seed = random_seed
  )

  # Combine all components into final dataset
  if (verbose) cat("- Combining all components...\n")
  final_data <- cbind(component1_data, component2_matrix, component3_matrix, component4_matrix)

  if (verbose) {
    cat(sprintf("- Dataset generated: %d observations, %d variables total\n",
                nrow(final_data), ncol(final_data) - 1))
    cat("- Component types:\n")
    cat("  1. GMMVar: Gaussian mixture with ascending class separation\n")
    cat("  2. NoiseGaussVar: Gaussian noise without class differences\n")
    cat("  3. UniformVar: Uniform with ascending class differences\n")
    cat("  4. NoiseUnifVar: Uniform noise regardless of class\n")
  }

  return(final_data)
}

#' Generate dataset with custom variable type proportions (original 4 distributions)
#'
#' @param random_seed Random seed for reproducibility
#' @param samples_per_class Number of samples per class
#' @param n_gmm Number of Gaussian mixture variables (Component 1)
#' @param n_gauss_noise Number of Gaussian noise variables (Component 2)
#' @param n_unif_noise Number of uniform noise variables (Component 3)
#' @param n_uniform Number of uniform variables with class differences (Component 4)
#' @param ... Additional parameters passed to component generation functions
#'
#' @return A data frame with specified proportions of variable types
generate_custom_synthetic_data <- function(random_seed = DEFAULT_RANDOM_SEED,
                                           samples_per_class = DEFAULT_SAMPLES_PER_CLASS,
                                           n_gmm = DEFAULT_NUM_VARIABLES,
                                           n_gauss_noise = DEFAULT_NUM_VARIABLES,
                                           n_uniform = DEFAULT_NUM_VARIABLES,
                                           n_unif_noise = DEFAULT_NUM_VARIABLES,
                                           ...) {

  # Component 1: Gaussian mixture data with ascending class separation
  component1_data <- generate_ascending_separation_data(
    n_samples = samples_per_class,
    n_vars = n_gmm,
    seed = random_seed,
    ...
  )

  # Component 2: Gaussian noise without class differences
  component2_matrix <- generate_gaussian_noise_data(
    n_samples = samples_per_class,
    n_vars = n_gauss_noise,
    seed = random_seed,
    ...
  )

  # Component 3: Uniform variables with ascending class differences
  component3_matrix <- generate_uniform_class_diff_data(
    n_samples = samples_per_class,
    n_vars = n_uniform,
    seed = random_seed,
    ...
  )

  # Component 4: Uniform random noise regardless of class
  component4_matrix <- generate_uniform_noise_data(
    n_samples = samples_per_class,
    n_vars = n_unif_noise,
    seed = random_seed,
    ...
  )

  # Combine all components
  final_data <- cbind(component1_data, component2_matrix, component3_matrix, component4_matrix)

  return(final_data)
}

#' Generate comprehensive biomedical synthetic dataset with all 10 distributions
#'
#' @param random_seed Random seed for reproducibility (default: DEFAULT_RANDOM_SEED)
#' @param samples_per_class Number of samples per class (default: DEFAULT_SAMPLES_PER_CLASS)
#' @param num_variables Number of variables for each distribution type (default: DEFAULT_NUM_VARIABLES)
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param ... Additional parameters passed to individual distribution functions
#'
#' @return A data frame with class labels and all 20 types of variables (10 with class differences + 10 noise)
generate_biomedical_synthetic_data <- function(random_seed = DEFAULT_RANDOM_SEED,
                                               samples_per_class = DEFAULT_SAMPLES_PER_CLASS,
                                               num_variables = DEFAULT_NUM_VARIABLES,
                                               verbose = TRUE,
                                               ...) {

  # Input validation
  if (samples_per_class <= 0 || num_variables <= 0) {
    stop("samples_per_class and num_variables must be positive integers")
  }

  if (verbose) {
    cat("Generating comprehensive biomedical synthetic dataset with 10 distributions...\n")
    cat(sprintf("- %d samples per class, %d variables per distribution type\n", samples_per_class, num_variables))
    cat("- Total variables: ", num_variables * 20, " (10 with class differences + 10 noise)\n")
  }

  # Initialize with class labels
  class_labels <- rep(1:2, each = samples_per_class)

  # Component list to store all matrices
  components <- list()

  # 1. Normal distribution
  if (verbose) cat("- Generating Normal distribution variables...\n")
  components$normal_diff <- generate_ascending_separation_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$normal_noise <- generate_gaussian_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 2. Uniform distribution
  if (verbose) cat("- Generating Uniform distribution variables...\n")
  components$uniform_diff <- generate_uniform_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$uniform_noise <- generate_uniform_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 3. Log-normal distribution
  if (verbose) cat("- Generating Log-normal distribution variables...\n")
  components$lognormal_diff <- generate_lognormal_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$lognormal_noise <- generate_lognormal_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 4. Binomial distribution
  if (verbose) cat("- Generating Binomial distribution variables...\n")
  components$binomial_diff <- generate_binomial_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$binomial_noise <- generate_binomial_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 5. Poisson distribution
  if (verbose) cat("- Generating Poisson distribution variables...\n")
  components$poisson_diff <- generate_poisson_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$poisson_noise <- generate_poisson_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 6. Exponential distribution
  if (verbose) cat("- Generating Exponential distribution variables...\n")
  components$exponential_diff <- generate_exponential_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$exponential_noise <- generate_exponential_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 7. Gamma distribution
  if (verbose) cat("- Generating Gamma distribution variables...\n")
  components$gamma_diff <- generate_gamma_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$gamma_noise <- generate_gamma_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 8. Beta distribution
  if (verbose) cat("- Generating Beta distribution variables...\n")
  components$beta_diff <- generate_beta_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$beta_noise <- generate_beta_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 9. Chi-squared distribution
  if (verbose) cat("- Generating Chi-squared distribution variables...\n")
  components$chisq_diff <- generate_chisq_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$chisq_noise <- generate_chisq_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # 10. Weibull distribution
  if (verbose) cat("- Generating Weibull distribution variables...\n")
  components$weibull_diff <- generate_weibull_class_diff_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)
  components$weibull_noise <- generate_weibull_noise_data(
    n_samples = samples_per_class, n_vars = num_variables, seed = random_seed, ...)

  # Combine all components into final dataset
  if (verbose) cat("- Combining all components...\n")

  # Start with class labels from the first component (normal_diff)
  final_data <- components$normal_diff

  # Add all other matrices (excluding the class column from matrices)
  for (i in 2:length(components)) {
    if (is.data.frame(components[[i]])) {
      # Remove class column if present
      var_cols <- setdiff(names(components[[i]]), "Cls")
      final_data <- cbind(final_data, components[[i]][, var_cols, drop = FALSE])
    } else {
      # It's a matrix
      final_data <- cbind(final_data, components[[i]])
    }
  }

  if (verbose) {
    cat(sprintf("- Dataset generated: %d observations, %d variables total\n",
                nrow(final_data), ncol(final_data) - 1))
    cat("- Distribution types (each with class differences and noise versions):\n")
    cat("  1. Normal distribution (GMMVar + NoiseGaussVar)\n")
    cat("  2. Uniform distribution (UniformVar + NoiseUnifVar)\n")
    cat("  3. Log-normal distribution (LogNormVar + NoiseLogNormVar)\n")
    cat("  4. Binomial distribution (BinomVar + NoiseBinomVar)\n")
    cat("  5. Poisson distribution (PoissonVar + NoisePoissonVar)\n")
    cat("  6. Exponential distribution (ExpVar + NoiseExpVar)\n")
    cat("  7. Gamma distribution (GammaVar + NoiseGammaVar)\n")
    cat("  8. Beta distribution (BetaVar + NoiseBetaVar)\n")
    cat("  9. Chi-squared distribution (ChiSqVar + NoiseChiSqVar)\n")
    cat("  10. Weibull distribution (WeibullVar + NoiseWeibullVar)\n")
  }

  return(final_data)
}

#' Generate fully customizable biomedical synthetic dataset
#'
#' @param random_seed Random seed for reproducibility
#' @param samples_per_class Number of samples per class
#' @param n_normal Number of normal variables (both class-diff and noise)
#' @param n_uniform Number of uniform variables (both class-diff and noise)
#' @param n_lognormal Number of log-normal variables (both class-diff and noise)
#' @param n_binomial Number of binomial variables (both class-diff and noise)
#' @param n_poisson Number of Poisson variables (both class-diff and noise)
#' @param n_exponential Number of exponential variables (both class-diff and noise)
#' @param n_gamma Number of gamma variables (both class-diff and noise)
#' @param n_beta Number of beta variables (both class-diff and noise)
#' @param n_chisq Number of chi-squared variables (both class-diff and noise)
#' @param n_weibull Number of Weibull variables (both class-diff and noise)
#' @param include_noise Whether to include noise versions of each distribution (default: TRUE)
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param ... Additional parameters passed to individual distribution functions
#'
#' @return A data frame with specified numbers of each distribution type
generate_custom_biomedical_data <- function(random_seed = DEFAULT_RANDOM_SEED,
                                            samples_per_class = DEFAULT_SAMPLES_PER_CLASS,
                                            n_normal = 0,
                                            n_uniform = 0,
                                            n_lognormal = 0,
                                            n_binomial = 0,
                                            n_poisson = 0,
                                            n_exponential = 0,
                                            n_gamma = 0,
                                            n_beta = 0,
                                            n_chisq = 0,
                                            n_weibull = 0,
                                            include_noise = TRUE,
                                            verbose = TRUE,
                                            ...) {

  # Input validation
  if (samples_per_class <= 0) {
    stop("samples_per_class must be positive integer")
  }

  # Count total variables
  total_vars <- n_normal + n_uniform + n_lognormal + n_binomial + n_poisson +
    n_exponential + n_gamma + n_beta + n_chisq + n_weibull

  if (include_noise) {
    total_vars <- total_vars * 2
  }

  if (total_vars == 0) {
    stop("At least one distribution type must have n > 0")
  }

  if (verbose) {
    cat("Generating custom biomedical synthetic dataset...\n")
    cat(sprintf("- %d samples per class, %d variables total\n", samples_per_class, total_vars))
  }

  # Initialize with class labels
  class_labels <- rep(1:2, each = samples_per_class)
  final_data <- data.frame(Cls = class_labels)

  # Generate each distribution type as requested
  if (n_normal > 0) {
    if (verbose) cat(sprintf("- Generating %d Normal variables...\n", n_normal * (1 + include_noise)))

    normal_diff <- generate_ascending_separation_data(
      n_samples = samples_per_class, n_vars = n_normal, seed = random_seed, ...)
    final_data <- cbind(final_data, normal_diff[, -1]) # Exclude Cls column

    if (include_noise) {
      normal_noise <- generate_gaussian_noise_data(
        n_samples = samples_per_class, n_vars = n_normal, seed = random_seed, ...)
      final_data <- cbind(final_data, normal_noise)
    }
  }

  if (n_uniform > 0) {
    if (verbose) cat(sprintf("- Generating %d Uniform variables...\n", n_uniform * (1 + include_noise)))

    uniform_diff <- generate_uniform_class_diff_data(
      n_samples = samples_per_class, n_vars = n_uniform, seed = random_seed, ...)
    final_data <- cbind(final_data, uniform_diff)

    if (include_noise) {
      uniform_noise <- generate_uniform_noise_data(
        n_samples = samples_per_class, n_vars = n_uniform, seed = random_seed, ...)
      final_data <- cbind(final_data, uniform_noise)
    }
  }

  if (n_lognormal > 0) {
    if (verbose) cat(sprintf("- Generating %d Log-normal variables...\n", n_lognormal * (1 + include_noise)))

    lognormal_diff <- generate_lognormal_class_diff_data(
      n_samples = samples_per_class, n_vars = n_lognormal, seed = random_seed, ...)
    final_data <- cbind(final_data, lognormal_diff)

    if (include_noise) {
      lognormal_noise <- generate_lognormal_noise_data(
        n_samples = samples_per_class, n_vars = n_lognormal, seed = random_seed, ...)
      final_data <- cbind(final_data, lognormal_noise)
    }
  }

  if (n_binomial > 0) {
    if (verbose) cat(sprintf("- Generating %d Binomial variables...\n", n_binomial * (1 + include_noise)))

    binomial_diff <- generate_binomial_class_diff_data(
      n_samples = samples_per_class, n_vars = n_binomial, seed = random_seed, ...)
    final_data <- cbind(final_data, binomial_diff)

    if (include_noise) {
      binomial_noise <- generate_binomial_noise_data(
        n_samples = samples_per_class, n_vars = n_binomial, seed = random_seed, ...)
      final_data <- cbind(final_data, binomial_noise)
    }
  }

  if (n_poisson > 0) {
    if (verbose) cat(sprintf("- Generating %d Poisson variables...\n", n_poisson * (1 + include_noise)))

    poisson_diff <- generate_poisson_class_diff_data(
      n_samples = samples_per_class, n_vars = n_poisson, seed = random_seed, ...)
    final_data <- cbind(final_data, poisson_diff)

    if (include_noise) {
      poisson_noise <- generate_poisson_noise_data(
        n_samples = samples_per_class, n_vars = n_poisson, seed = random_seed, ...)
      final_data <- cbind(final_data, poisson_noise)
    }
  }

  if (n_exponential > 0) {
    if (verbose) cat(sprintf("- Generating %d Exponential variables...\n", n_exponential * (1 + include_noise)))

    exponential_diff <- generate_exponential_class_diff_data(
      n_samples = samples_per_class, n_vars = n_exponential, seed = random_seed, ...)
    final_data <- cbind(final_data, exponential_diff)

    if (include_noise) {
      exponential_noise <- generate_exponential_noise_data(
        n_samples = samples_per_class, n_vars = n_exponential, seed = random_seed, ...)
      final_data <- cbind(final_data, exponential_noise)
    }
  }

  if (n_gamma > 0) {
    if (verbose) cat(sprintf("- Generating %d Gamma variables...\n", n_gamma * (1 + include_noise)))

    gamma_diff <- generate_gamma_class_diff_data(
      n_samples = samples_per_class, n_vars = n_gamma, seed = random_seed, ...)
    final_data <- cbind(final_data, gamma_diff)

    if (include_noise) {
      gamma_noise <- generate_gamma_noise_data(
        n_samples = samples_per_class, n_vars = n_gamma, seed = random_seed, ...)
      final_data <- cbind(final_data, gamma_noise)
    }
  }

  if (n_beta > 0) {
    if (verbose) cat(sprintf("- Generating %d Beta variables...\n", n_beta * (1 + include_noise)))

    beta_diff <- generate_beta_class_diff_data(
      n_samples = samples_per_class, n_vars = n_beta, seed = random_seed, ...)
    final_data <- cbind(final_data, beta_diff)

    if (include_noise) {
      beta_noise <- generate_beta_noise_data(
        n_samples = samples_per_class, n_vars = n_beta, seed = random_seed, ...)
      final_data <- cbind(final_data, beta_noise)
    }
  }

  if (n_chisq > 0) {
    if (verbose) cat(sprintf("- Generating %d Chi-squared variables...\n", n_chisq * (1 + include_noise)))

    chisq_diff <- generate_chisq_class_diff_data(
      n_samples = samples_per_class, n_vars = n_chisq, seed = random_seed, ...)
    final_data <- cbind(final_data, chisq_diff)

    if (include_noise) {
      chisq_noise <- generate_chisq_noise_data(
        n_samples = samples_per_class, n_vars = n_chisq, seed = random_seed, ...)
      final_data <- cbind(final_data, chisq_noise)
    }
  }

  if (n_weibull > 0) {
    if (verbose) cat(sprintf("- Generating %d Weibull variables...\n", n_weibull * (1 + include_noise)))

    weibull_diff <- generate_weibull_class_diff_data(
      n_samples = samples_per_class, n_vars = n_weibull, seed = random_seed, ...)
    final_data <- cbind(final_data, weibull_diff)

    if (include_noise) {
      weibull_noise <- generate_weibull_noise_data(
        n_samples = samples_per_class, n_vars = n_weibull, seed = random_seed, ...)
      final_data <- cbind(final_data, weibull_noise)
    }
  }

  if (verbose) {
    cat(sprintf("- Dataset generated: %d observations, %d variables total\n",
                nrow(final_data), ncol(final_data) - 1))
  }

  return(final_data)
}