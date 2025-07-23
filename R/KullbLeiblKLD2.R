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
  num_unique_x <- length(x)

  # If there are less than 2 unique values, return 0 for the KLD
  if (num_unique_x < 2) {
    return(list(KLD = 0, p = 1, q = 1, x = x))
  }

  # Calculate the frequency distributions
  freq_p <- tabulate(match(P, x))
  freq_q <- tabulate(match(Q, x))
  p <- freq_p / sum(freq_p)
  q <- freq_q / sum(freq_q)

  # Calculate the Kullback-Leibler divergence
  EPS <- 1 / 10000
  log_p <- log(p[p > EPS])
  log_q <- log(q[q > EPS])
  log_p_to_q <- log(p[p > EPS] / q[q > EPS])
  log_q_to_p <- log(q[q > EPS] / p[p > EPS])

  KLDpq <- sum(log_p * log_p_to_q, na.rm = TRUE) / num_unique_x
  KLDqp <- sum(log_q * log_q_to_p, na.rm = TRUE) / num_unique_x

  if (sym) {
    KLD <- KLDqp + KLDpq
  } else {
    KLD <- KLDpq
  }

  return(list(KLD = KLD, p = p, q = q, x = x))
}