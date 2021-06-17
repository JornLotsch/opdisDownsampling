\name{opdisDownsampling}
\alias{opdisDownsampling}
\title{Optimal Distribution Preserving Down-Sampling of Bio-Medical Data}
\description{
  The package provides the necessary functions for optimal distribution-preserving down-sampling of large (bio-medical) data sets.
}
\usage{
opdisDownsampling(Data, Cls, Size, Seed, nTrials,
  TestStat = "ad", MaxCores = 8, JobSize = 10000, PCAimportance = FALSE)
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
 \item{PCAimportance}{PCA based feature selection; only variables important in PCA projection are considered.}
}
\value{
  Returns a list of data containing the drawn samples and the omitted data.
 \item{ReducedData}{the selected sample data and class information.}
 \item{RemucedData}{the not-selected sample data and class information.}
 }
\references{
 Lotsch, J. (2021):
 opdisDownsampling -- an R Package for optimal distribution-preserving down-sampling of large (bio-medical) data sets.
}
\author{
  Jorn Lotsch
}
\examples{
## example 1
data(iris)
Iris50percent <- opdisDownsampling(Data = iris[,1:4], Cls = as.integer(iris$Species),
  Size = 50, MaxCores = 1)
}
\keyword{opdisDownsampling}
\keyword{data sampling}
