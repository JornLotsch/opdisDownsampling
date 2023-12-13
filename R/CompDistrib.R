# Several methods to compare two data distributions
#' @importFrom twosamples ad_stat kuiper_stat cvm_stat wass_stat dts_stat
#' @importFrom stats ks.test prcomp na.omit
#' @importFrom EucDist EucDist
#' @importFrom KullbLeiblKLD2 KullbLeiblKLD2
#' @importFrom amrdd amrdd
CompDistrib <- function(vector1, vector2, TestStat) {
  if (length(vector1[!is.na(vector1)]) == 0 |
      length(vector2[!is.na(vector2)]) == 0) {
    Stati <- 1e+27
  } else {
    Stati <- switch(
      TestStat,
      ad = {
        twosamples::ad_stat(vector1, vector2)
      },
      kuiper = {
        twosamples::kuiper_stat(na.omit(vector1), na.omit(vector2))
      },
      cvm = {
        twosamples::cvm_stat(na.omit(vector1), na.omit(vector2))
      },
      wass = {
        twosamples::wass_stat(na.omit(vector1), na.omit(vector2))
      },
      dts = {
        twosamples::dts_stat(na.omit(vector1), na.omit(vector2))
      },
      ks = {
        ks.test(na.omit(vector1), na.omit(vector2))$statistic
      },
      kld = {
        KullbLeiblKLD2(na.omit(vector1), na.omit(vector2))$KLD
      },
      amrdd = {
        amrdd(na.omit(vector1), na.omit(vector2))
      },
      euc = {
        EucDist(na.omit(vector1), na.omit(vector2))
      }
    )
  }
  return(Stati)
}
