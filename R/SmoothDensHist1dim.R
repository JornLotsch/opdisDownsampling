#Function to calculate a smooth density histogram for a single variable
#' @importFrom pracma histc linspace
#' @importFrom caTools trapz

nanmax <- function (Data) {
  if (length(dim(Data)) == 2) {
    SpaltenMinima <- apply(Data, 2, function(x) max(x, na.rm = TRUE))
    SpaltenInd <- NaN
  } else {
    SpaltenMinima <- max(Data, na.rm = TRUE)
    SpaltenInd <- which(Data == SpaltenMinima)
  }
  return(SpaltenMinima)
}

SmoothDensHist1dim <- function (x, KernelsOrNbins = NULL, SDHinPercent, lambda) {
  if (length(x) == 0) {
    warning("SmoothDensHist1dim: Size of x is zero.", call. = FALSE)
    if (is.null(KernelsOrNbins) == TRUE) Kernels = 1 else Kernels = KernelsOrNbins
    SDH = Kernels * 0
  } else {
    requireNamespace("pracma")
    requireNamespace("caTools")

    smooth1D = function(Y, lambda) {
      if (is.vector(Y))
        Y = as.matrix(Y)
      Y[is.na(Y)] <-0
      dd = dim(Y)
      m = dd[1]
      n = dd[2]
      E = diag(m)
      D1 = (diff(E, 1))
      D2 = (diff(D1, 1))
      P = lambda ^ 2 * (t(D2) %*% D2) + 2 * lambda * (t(D1) %*% D1)
      Z = solve((E + P), Y)
      return(Z)
    }
    if (missing(lambda)) lambda = 20
    if (missing(SDHinPercent)) SDHinPercent = 0
    if (is.null(KernelsOrNbins))  KernelsOrNbins = 200
    if (length(KernelsOrNbins) < 1) KernelsOrNbins = 200

    x = x[is.finite(x)]
    n = length(x)
    minx = min(x, na.rm = T)
    maxx = max(x, na.rm = T)

    if (length(KernelsOrNbins) == 1) {
      nbins = KernelsOrNbins
      edges1 = seq(from = minx, to = maxx, length.out = (nbins + 1))
      end = length(edges1)
      Kernels = c(edges1[1:(end - 1)] + 0.5 * diff(edges1))
      InInd = c()
    } else {
      Kernels = c(KernelsOrNbins)
    }

    InInd = which((Kernels >= minx) & (Kernels <= maxx))

    if (length(InInd) == 0) {
      SDH = Kernels * 0
      nbins = 1
    } else {
      DataInd = which((x >= min(Kernels[InInd])) & (x <= max(Kernels[InInd])))
      if (length(DataInd) < 2) {
        if (length(DataInd) == 0) {
          SDH = Kernels * 0
          nbins = 1
        } else {
          SDH = Kernels * 0
          SDH[InInd] <- 1
          nbins = 1
        }
      } else {
        x = x[DataInd]
        edges1 = Kernels[InInd]
        nbins = length(edges1)
        #Kernels <- edges1
        end = length(edges1)
        edges1 = c(-Inf, edges1[2:(end - 1)], Inf)
        V = pracma::histc(x, edges1)
        dummy = V$cnt
        bin = V$bin
        H = dummy / n
        SDH = smooth1D(H, nbins / lambda)
        SDH = as.vector(SDH)
      }

      if (length(DataInd) > 1) {
        sdh = SDH
        SDH = Kernels * 0
        SDH[InInd] = sdh
      }
    }

    if (SDHinPercent == 0) {
      if (sum(SDH) == 0) {
        Area = 0
      } else {
        Area = caTools::trapz(Kernels, SDH)
      }
      if (Area < 1e-10) {
        SDH = rep(0, length(Kernels))
      } else {
        SDH = SDH / Area
      }
    } else {
      SDH = SDH / nanmax(SDH)
    }
  }
  return(list(Kernels = Kernels, SDH = SDH))
}
