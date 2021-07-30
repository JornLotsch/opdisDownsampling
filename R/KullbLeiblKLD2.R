# Function to calculate the Kullback-leibler divergence
KullbLeiblKLD2 <- function(P, Q, sym = 1) {
  EPS <- 1/10000
  P <- P[!is.na(P)]
  Q <- Q[!is.na(Q)]
  x <- sort(na.last = T, unique(c(P, Q)))
  AnzUniqX <- length(x)
  if (AnzUniqX < 2) {
    KLD <- 0
    KLDpq <- 0
    KLDqp <- 0
    p <- 1
    q <- 1
  } else {
    FreqP <- vector()
    for (i in 1:length(x)) FreqP[i] <- length(which(P == x[i]))
    FreqQ <- vector()
    for (i in 1:length(x)) FreqQ[i] <- length(which(Q == x[i]))
    p <- FreqP/sum(FreqP)
    q <- FreqQ/sum(FreqQ)
    LogP <- vector("numeric", AnzUniqX) * NaN
    LogQ <- vector("numeric", AnzUniqX) * NaN
    LogPzuQ <- vector("numeric", AnzUniqX) * NaN
    LogQzuP <- vector("numeric", AnzUniqX) * NaN
    Ind <- which(p > EPS, arr.ind = T)
    LogP[Ind] <- log(p[Ind])
    Ind <- which(q > EPS, arr.ind = T)
    LogQ[Ind] <- log(q[Ind])
    Ind <- which((q > EPS) & (p > EPS), arr.ind = T)
    LogPzuQ[Ind] <- log(p[Ind]/q[Ind])
    Ind <- which((q > EPS) & (p > EPS), arr.ind = T)
    LogQzuP[Ind] <- log(q[Ind]/p[Ind])
    KLDpq <- LogP * LogPzuQ
    ZeroInd <- which(p <= EPS, arr.ind = T)
    KLDpq[ZeroInd] <- 0
    KLDqp <- LogQ * LogQzuP
    ZeroInd <- which(q <= EPS, arr.ind = T)
    KLDqp[ZeroInd] <- 0
    KLDpq[!is.finite(KLDpq)] <- NaN
    KLDqp[!is.finite(KLDqp)] <- NaN
    KLDpqS <- sum(KLDpq[is.finite(KLDpq)], na.rm = T)/AnzUniqX
    KLDqpS <- sum(KLDqp[is.finite(KLDqp)], na.rm = T)/AnzUniqX
    if (sym == 1) {
      KLD <- KLDqpS + KLDpqS
    } else {
      KLD <- KLDpqS
    }
  }
  return(list(KLD <- KLD, p <- p, q <- q, x <- x))
}
