\name{opdisDownsampling}
\alias{opdisDownsampling}
\title{Distribution preserving downsampling of biomedical data}
\description{
  The package provides the necessary functions for class proportinal daowsampling data.
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
  Returns a list of data including the drwan samples and the left-out data.
 \item{ReducedData}{the selected sample data and class information.}
 \item{RemucedData}{the not selected data and class information.}
 }
\references{
 Lotsch, J. (2021): 
 opdisDownsampling -- an R Package for distribution preserving downsampling of biomedical data.
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
