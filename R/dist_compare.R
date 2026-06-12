#' Compare Two Data Distributions
#'
#' This internal helper compares two data distributions using one of the
#' supported distribution-similarity statistics.
#'
#' @param vector1 A numeric vector representing the first data distribution.
#' @param vector2 A numeric vector representing the second data distribution.
#' @param TestStat A character string specifying the statistic to use. Supported
#'   values are \code{"ad"}, \code{"kuiper"}, \code{"cvm"}, \code{"wass"},
#'   \code{"dts"}, \code{"ks"}, \code{"kld"}, \code{"amrdd"}, \code{"euc"},
#'   and \code{"nent"}.
#'
#' @return The value of the specified distribution statistic.
#'
#' @importFrom twosamples ad_stat kuiper_stat cvm_stat wass_stat dts_stat
#' @importFrom stats ks.test na.omit
#'
CompDistrib <- function(vector1, vector2, TestStat) {
  # Check if both input vectors have at least one non-NA value
  if (length(vector1[!is.na(vector1)]) == 0 || length(vector2[!is.na(vector2)]) == 0) {
    return(1e+27)
  }

  # Remove NA values from the input vectors
  vector1 <- na.omit(vector1)
  vector2 <- na.omit(vector2)

  # Perform the specified statistical test
  Stati <- switch(TestStat,
    ad = twosamples::ad_stat(vector1, vector2),
    kuiper = twosamples::kuiper_stat(vector1, vector2),
    cvm = twosamples::cvm_stat(vector1, vector2),
    wass = twosamples::wass_stat(vector1, vector2),
    dts = twosamples::dts_stat(vector1, vector2),
    ks = ks.test(vector1, vector2)$statistic,
    kld = KullbLeiblKLD2(vector1, vector2)$KLD,
    amrdd = amrdd(vector1, vector2),
    euc = EucDist(vector1, vector2),
    nent = abs_norm_entropy_diff(vector1, vector2),

    # If the specified test is not supported, return a large value
    stop("Unsupported test statistic: ", TestStat)
  )

  return(Stati)
}
