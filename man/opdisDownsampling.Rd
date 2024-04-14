\name{opdisDownsampling}
\alias{opdisDownsampling}
\title{Optimal Distribution Preserving Down-Sampling of Bio-Medical Data}
\description{
  The package provides the necessary functions for optimal distribution-preserving down-sampling of large (bio-medical) data sets.
}
\usage{
opdisDownsampling(Data, Cls, Size, Seed, nTrials = 1000,
TestStat = "ad", MaxCores = getOption("mc.cores", 2L), PCAimportance = FALSE)
}
\arguments{
 \item{Data}{the (numerical!) data as a vector, matrix or data frame.}
 \item{Cls}{the class information, if any, as a vector of similar length as instances in the data.}
 \item{Size}{the total number of instances across all classes to be drawn.}
 \item{Seed}{a predefined seed to modify the results.}
 \item{nTrials}{how many samples to choose from should be randomly drawn.}
 \item{TestStat}{statistical criterion for similarity judgment.}
 \item{MaxCores}{maximum number of cpu cores to use for parallel computing.}
 \item{PCAimportance}{PCA based feature selection; only variables important in PCA projection are considered.}
}
\value{
  Returns a list of data containing the drawn samples and the omitted data.
 \item{ReducedData}{the selected sample data and class information.}
 \item{ReducedData}{the not-selected sample data and class information.}
 \item{ReducedInstances}{the instance numbers of the selected sample data.}
 }
\references{
 Lotsch, J., Malkusch, S., Ultsch, A. (2021):
 Optimal distribution-preserving downsampling of large biomedical data sets (opdisDownsampling).
 PLoS One. 2021 Aug 5;16(8):e0255838. doi: 10.1371/journal.pone.0255838. eCollection 2021.
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
