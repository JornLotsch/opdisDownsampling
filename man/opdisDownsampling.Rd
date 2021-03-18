\name{opdisDownsampling}
\alias{opdisDownsampling}
\title{Optimal distribution-preserving downsampling of (biomedical) data}
\description{
  The package provides the necessary functions for optimal distribution-preserving downsampling of large (biomedical) data sets.
}
\usage{
opdisDownsampling(Data, Cls, Size, Seed, nTrials,
  TestStat = "ad", MaxCores = 8, JobSize = 10000)
}
\arguments{
 \item{Data}{the data as a vector, matrix or data frame.}
 \item{Cls}{the class information, if any, as a vector of similar length as instances in the data.}
 \item{Size}{the percentage of instances per class to be drwan.}
 \item{Seed}{a predefined seed to modify the results.}
 \item{nTrials}{how many samples to choose from should be randomly drawn.}
 \item{TestStat}{statistical criterion for similarity judgment.}
 \item{MaxCores}{maximum number of cpu cores to use for parallel computing.}
 \item{JobSize}{how many samples can be drwan at once.}
}
\value{
  Returns a list of data containing the drawn samples and the omitted data.
 \item{ReducedData}{the selected sample data and class information.}
 \item{RemucedData}{the not-selected sample data and class information.}
 }
\references{
 Lotsch, J. (2021):
 opdisDownsampling -- an R Package for optimal distribution-preserving downsampling of large (biomedical) data sets.
}
\author{
  Jorn Lotsch
}
\examples{
## example 1
data(iris)
Iris50percent <- opdisDownsampling(iris[,1:4], as.integer(iris$Species), 50)
}
\keyword{opdisDownsampling}
\keyword{data sampling}
