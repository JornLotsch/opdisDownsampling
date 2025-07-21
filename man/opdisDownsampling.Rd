\name{opdisDownsampling}
\alias{opdisDownsampling}
\title{Optimal Distribution Preserving Down-Sampling of Bio-Medical Data}
\description{
  The package provides functions for optimal distribution-preserving down-sampling of large (bio-medical) data sets. It draws statistically representative subsets of data while preserving the class proportions and original data distribution.
}
\usage{
opdisDownsampling(Data, Cls, Size, Seed, nTrials = 1000,
  TestStat = "ad", MaxCores = getOption("mc.cores", 2L),
  PCAimportance = FALSE, CheckRemoved = FALSE,
  JobSize = NULL, verbose = FALSE)
}
\arguments{
 \item{Data}{Numeric data as a vector, matrix, or data frame. Each row represents an instance, each column a variable.}

 \item{Cls}{Optional vector with class labels for each instance in \code{Data}. If missing, all instances are treated as belonging to the same class.}

 \item{Size}{The number (integer) or proportion (0<Size<1) of instances to draw from the dataset.
      The reduction is class proportional and aims to preserve the variable distributions.}

 \item{Seed}{Integer seed used to initialize the random number generator. Useful for reproducibility.}

 \item{nTrials}{Number of random sampling trials used to find the optimal subset (default: 1000).}

 \item{TestStat}{Character string defining the statistical test used to assess distribution similarity.

 Available options are:
   \itemize{
     \item{\code{"ad"}: Anderson–Darling statistic}
     \item{\code{"kuiper"}: Kuiper statistic}
     \item{\code{"cvm"}: Cramér–von Mises statistic}
     \item{\code{"wass"}: Wasserstein distance}
     \item{\code{"dts"}: Distributional Transform Statistic}
     \item{\code{"ks"}: Kolmogorov–Smirnov statistic}
     \item{\code{"kld"}: Kullback–Leibler divergence (via \code{KullbLeiblKLD2()})}
     \item{\code{"amrdd"}: Average Mean Root of Distributional Differences (via \code{amrdd()})}
     \item{\code{"euc"}: Euclidean distance (via \code{EucDist()})}
   }
 }

 \item{MaxCores}{Maximum number of CPU cores to use for parallel computing (default is value stored in \code{getOption("mc.cores")}, or 2 if missing).}

 \item{PCAimportance}{Logical; if \code{TRUE}, only variables deemed important by principal component analysis are used in computing similarity statistics.}

 \item{CheckRemoved}{Logical; if \code{TRUE}, the optimization also includes the subset of removed cases,
      which can be used for data splits into training/test and validation.}

 \item{JobSize}{Integer specifying the number of trials to process in each memory-efficient chunk.
      If \code{NULL} (default), an optimal chunk size is automatically calculated based on data dimensions,
      number of trials, available system memory, and number of processor cores. Smaller chunk sizes use less
      memory but may increase computational overhead, while larger chunk sizes are more efficient but require
      more memory. Manual specification overrides the automatic calculation.}

 \item{verbose}{Logical; if \code{TRUE}, prints diagnostic information about chunk size calculation,
      including data dimensions, estimated memory usage, and chosen chunking strategy. Useful for
      understanding memory usage patterns and debugging performance issues (default: \code{FALSE}).}
}
\value{
  Returns a list with the following elements:
  \item{ReducedData}{Down-sampled data set (as data frame or matrix) including only the selected instances.}
  \item{RemovedData}{Data not included in the sample.}
  \item{ReducedInstances}{Row indices (or names) of the selected instances from the original data set.}
  \item{RemovedInstances}{Row indices (or names) of the unselected instances from the original data set.}
}
\details{
  The function automatically optimizes memory usage through chunked processing when dealing with large datasets
  or high numbers of trials. The chunking strategy considers:
  \itemize{
    \item{Data size (number of rows and columns)}
    \item{Available system memory (detected on Linux systems)}
    \item{Number of processor cores}
    \item{Number of trials to perform}
  }

  For very large datasets or memory-constrained systems, smaller chunks are automatically selected to prevent
  memory exhaustion, while smaller datasets use larger chunks to minimize computational overhead.
}
\references{
  Lotsch, J., Malkusch, S., Ultsch, A. (2021):\\
  Optimal distribution-preserving downsampling of large biomedical data sets.\\
  \emph{PLoS ONE} 16(8): e0255838. \doi{10.1371/journal.pone.0255838}
}
\author{
  Jorn Lotsch
}
\examples{
## Example: Down-sample the Iris dataset to 50 points
data(iris)
Iris50percent <- opdisDownsampling(Data = iris[,1:4], Cls = as.integer(iris$Species),
  Size = 50, MaxCores = 1)

## Example: Down-sample with custom chunk size and verbose output
data(iris)
Iris50percent <- opdisDownsampling(Data = iris[,1:4], Cls = as.integer(iris$Species),
  Size = 50, MaxCores = 1, JobSize = 25, verbose = TRUE)

## Example: Memory-efficient processing of large dataset with many trials
\dontrun{
# For large datasets, automatic chunking optimizes memory usage
LargeDataSample <- opdisDownsampling(Data = large_dataset,
  Size = 0.1, nTrials = 5000, verbose = TRUE)
}
}
\keyword{downsampling}
\keyword{data reduction}
\keyword{data splits}
\keyword{distribution preservation}
\keyword{sampling}
\keyword{memory optimization}