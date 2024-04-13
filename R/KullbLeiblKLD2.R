#' Calculate the Kullback-Leibler Divergence
#'
#' This function calculates the Kullback-Leibler divergence between two probability
#' distributions.
#'
#' @param P A numeric vector representing the first probability distribution.
#' @param Q A numeric vector representing the second probability distribution.
#' @param sym A logical value indicating whether to use the symmetric Kullback-Leibler
#'   divergence (TRUE) or the standard Kullback-Leibler divergence (FALSE).
#'
#' @return A list with the following elements:
#'   - `KLD`: The Kullback-Leibler divergence between the two distributions.
#'   - `p`: The normalized frequency distribution of the first distribution.
#'   - `q`: The normalized frequency distribution of the second distribution.
#'   - `x`: The unique values in the combined distribution.
#'
KullbLeiblKLD2 <- function(P, Q, sym = TRUE) {
  # Remove NA values from the input vectors
  P <- P[!is.na(P)]
  Q <- Q[!is.na(Q)]
  
  # Combine the unique values from both distributions
  x <- sort(unique(c(P, Q)), na.last = TRUE)
  AnzUniqX <- length(x)
  
  # If there are less than 2 unique values, return 0 for the KLD
  if (AnzUniqX < 2) {
    return(list(KLD = 0, p = 1, q = 1, x = x))
  }
  
  # Calculate the frequency distributions
  FreqP <- tabulate(match(P, x))
  FreqQ <- tabulate(match(Q, x))
  p <- FreqP / sum(FreqP)
  q <- FreqQ / sum(FreqQ)
  
  # Calculate the Kullback-Leibler divergence
  EPS <- 1 / 10000
  LogP <- log(p[p > EPS])
  LogQ <- log(q[q > EPS])
  LogPzuQ <- log(p[p > EPS] / q[q > EPS])
  LogQzuP <- log(q[q > EPS] / p[p > EPS])
  
  KLDpq <- sum(LogP * LogPzuQ, na.rm = TRUE) / AnzUniqX
  KLDqp <- sum(LogQ * LogQzuP, na.rm = TRUE) / AnzUniqX
  
  if (sym) {
    KLD <- KLDqp + KLDpq
  } else {
    KLD <- KLDpq
  }
  
  return(list(KLD = KLD, p = p, q = q, x = x))
}
