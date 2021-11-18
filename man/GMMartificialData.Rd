\name{GMMartificialData}
\alias{GMMartificialData}

\docType{data}
\title{Example data an artificial Gaussian mixture.}

\description{Dataset of 30000 instances with 10 variables that are Gaussian mixtures and belong to classes Cls = 1, 2, or 3, with different means and standard deviations and equal weights of 0.5, 0.4, and 0.1, respectively.}

\usage{data("GMMartificialData")}

\details{
Size 30000 x 10, stored in \code{GMMartificialData$[X1,X2,X3,X4,X5,X6,X7,X8,X9,X10]}

Classes 3, stored in \code{GMMartificialData$Cls}
}

\examples{
data(GMMartificialData)
str(GMMartificialData)
}
