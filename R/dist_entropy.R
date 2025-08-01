#' Absolute Normalized Entropy Difference Between Two Data Vectors
#'
#' This function computes the absolute difference in normalized Shannon entropies
#' between two input vectors (numeric, factor, or character), using adaptive
#' binning for numeric types and standard normalization by maximal entropy.
#'
#' @param vector1 A numeric, integer, factor, or character vector.
#' @param vector2 A numeric, integer, factor, or character vector.
#' @param bins_method A binning method for numeric vectors (default "FD" for Freedman-Diaconis; see ?hist).
#' @param min_bins Minimum number of bins for numeric binning (default = 2).
#' @param na_as_category Logical: should NA be considered as a special category? (default = TRUE).
#'
#' @return The absolute difference in normalized entropy between the two vectors (range: 0 to 1).
#'
#' @examples
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- runif(100)
#' abs_norm_entropy_diff(x1, x2)
#'
abs_norm_entropy_diff <- function(vector1, vector2, bins_method = "FD", min_bins = 2, na_as_category = TRUE) {

  normalized_entropy <- function(x, bins_method = "FD", min_bins = 2, na_as_category = TRUE) {
    if (is.numeric(x) || is.integer(x)) {
      x_no_na <- na.omit(x)
      # Adaptive binning
      bins <- tryCatch(
        hist(x_no_na, breaks = bins_method, plot = FALSE)$breaks,
        error = function(e) seq(min(x_no_na), max(x_no_na), length.out = min_bins + 1)
      )
      if (length(bins) - 1 < min_bins) {
        bins <- seq(min(x_no_na), max(x_no_na), length.out = min_bins + 1)
      }
      x_disc <- cut(x, breaks = bins, include.lowest = TRUE, right = FALSE, labels = FALSE)
      if (na_as_category) x_disc[is.na(x_disc)] <- max(x_disc, na.rm = TRUE) + 1
      p <- prop.table(table(x_disc))
    } else {
      if (na_as_category) x[is.na(x)] <- "NA"
      p <- prop.table(table(x))
    }
    H <- -sum(p * log(p))
    n <- length(p)
    H_max <- if (n > 1) log(n) else 1 # Avoid division by zero
    H_norm <- H / H_max
    return(H_norm)
  }

  h1 <- normalized_entropy(vector1, bins_method, min_bins, na_as_category)
  h2 <- normalized_entropy(vector2, bins_method, min_bins, na_as_category)
  abs(h1 - h2)
}
