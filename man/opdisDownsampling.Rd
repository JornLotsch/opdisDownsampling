\name{opdisDownsampling}
\alias{opdisDownsampling}
\title{Optimal Distribution Preserving Down-Sampling of Bio-Medical Data}
\description{
  The package provides functions for optimal distribution-preserving down-sampling of large (bio-medical) data sets.
  It draws statistically representative subsets of data while preserving the class proportions and original data distribution.
}
\usage{
opdisDownsampling(Data, Cls, Size, Seed = "simple", nTrials = 1000,
  TestStat = "ad", MaxCores = getOption("mc.cores", 2L),
  PCAimportance = FALSE, CheckRemoved = FALSE, CheckThreefold = FALSE,
  OptimizeBetween = FALSE, JobSize = 0, verbose = FALSE,
  NonNoiseSelection = FALSE, UniformTestStat = "ks",
  UniformThreshold = 0.05, WorstSample = FALSE)
}
\arguments{
 \item{Data}{Numeric data as a vector, matrix, or data frame. Each row represents an instance, each column a variable.}

 \item{Cls}{Optional vector with class labels for each instance in \code{Data}. If missing, all instances are treated as belonging to the same class.}

 \item{Size}{The number (integer) or proportion (0<Size<1) of instances to draw from the dataset.
      The reduction is class proportional and aims to preserve the variable distributions.}

 \item{Seed}{Seed value for reproducibility. Options: \code{"auto"} (complex seed recovery, slow), \code{"simple"}
      (fast reproducible seed, default), or integer (exact control, recommended). Use integers for systematic testing with different seeds,
      \code{"simple"} for general use, or \code{"auto"} when you need to maintain exact RNG state continuity.}

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
     \item{\code{"nent"}: Absolute normalized entropy difference (via \code{abs_norm_entropy_diff()})}   
     }
 }

 \item{MaxCores}{Maximum number of CPU cores to use for parallel computing (default is value stored in \code{getOption("mc.cores")}, or 2 if missing).}

 \item{PCAimportance}{Logical; if \code{TRUE}, only variables deemed important by principal component analysis are used in computing similarity statistics.}

 \item{CheckRemoved}{Logical; if \code{TRUE}, the optimization also includes the subset of removed cases,
      which can be used for data splits into training/test and validation.}

 \item{CheckThreefold}{Logical; if \code{TRUE}, the optimization also includes the equal distribution of
      the reduced and the removed data subsets, e.g. training/test and validation. Ignored when CheckRemoved is FALSE}

 \item{OptimizeBetween}{Logical; if \code{TRUE}, the optimization targets the reduced versus the removed data subsets, e.g. training/test and validation.
      If set, all other comparisons are not performed.}

 \item{JobSize}{Integer specifying the number of trials to process in each memory-efficient chunk.
      If \code{0} (default), no chunking is applied. If \code{NULL}, an optimal chunk size is automatically calculated based on data dimensions,
      number of trials, available system memory, and number of processor cores. Smaller chunk sizes use less
      memory but may increase computational overhead, while larger chunk sizes are more efficient but require
      more memory. Manual specification overrides the automatic calculation.}

 \item{verbose}{Logical; if \code{TRUE}, prints diagnostic information about chunk size calculation,
      including data dimensions, estimated memory usage, and chosen chunking strategy. Useful for
      understanding memory usage patterns and debugging performance issues (default: \code{FALSE}).}

 \item{NonNoiseSelection}{Logical; if \code{TRUE}, uses statistical distribution tests to identify
      and select only variables that significantly deviate from uniform distributions. This feature
      helps filter out noise variables and focus optimization on meaningful signal variables
      (default: \code{FALSE}).}

 \item{UniformTestStat}{Character string specifying the statistical test for non-uniform variable
      selection. Uses the same test statistics as \code{TestStat}. Only used when
      \code{NonNoiseSelection = TRUE} (default: \code{"ks"}).}

 \item{UniformThreshold}{Numeric threshold value above which variables are considered significantly
      non-uniform. Higher values make the selection more stringent. Only used when
      \code{NonNoiseSelection = TRUE} (default: \code{0.1}).}

 \item{WorstSample}{Logical; if \code{TRUE}, the optimization targets the least similar subsample.
       For testing  purpose (default: FALSE).}
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

  \strong{Variable Selection Methods:}

  Two variable selection methods can be used individually or in combination:
  \itemize{
    \item{\strong{PCA-based selection} (\code{PCAimportance = TRUE}): Identifies variables with high loadings
          in the first principal components, focusing on variables that capture the most variance in the data.}
    \item{\strong{Non-uniform selection} (\code{NonNoiseSelection = TRUE}): Identifies variables that
          significantly deviate from uniform distributions, helping to filter out noise variables and
          focus on meaningful signal variables.}
  }

  When both methods are enabled, the function first attempts to use their intersection. If no variables
  are found in the intersection, it uses their union. If no variables are selected by any method,
  the first variable is used to ensure the analysis can proceed, with a warning message displayed.
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
  Size = 50, Seed = 42, MaxCores = 1)

## Example: Down-sample with custom chunk size and verbose output
data(iris)
Iris50percent <- opdisDownsampling(Data = iris[,1:4], Cls = as.integer(iris$Species),
  Size = 50, Seed = 42, MaxCores = 1, JobSize = 25, verbose = TRUE)

## Example: Use non-uniform variable selection to filter noise variables
data(iris)
# Add some uniform noise variables to demonstrate filtering
noisy_iris <- cbind(iris[,1:4],
                    noise1 = runif(nrow(iris), 0, 10),
                    noise2 = runif(nrow(iris), -5, 5))
Iris_filtered <- opdisDownsampling(Data = noisy_iris, Cls = as.integer(iris$Species),
  Size = 50, Seed = 42, NonNoiseSelection = TRUE, UniformThreshold = 0.05, MaxCores = 1)

## Example: Combine PCA and non-uniform selection
data(iris)
Iris_combined <- opdisDownsampling(Data = iris[,1:4], Cls = as.integer(iris$Species),
  Size = 50, Seed = 42, PCAimportance = TRUE, NonNoiseSelection = TRUE, MaxCores = 1)

## Example: Memory-efficient processing of large dataset with many trials
\dontrun{
# For large datasets, automatic chunking optimizes memory usage
LargeDataSample <- opdisDownsampling(Data = large_dataset,
  Size = 0.1, Seed = 42, nTrials = 5000, verbose = TRUE)
}
}
\keyword{downsampling}
\keyword{data reduction}
\keyword{data splits}
\keyword{distribution preservation}
\keyword{sampling}
\keyword{memory optimization}
\keyword{variable selection}
\keyword{noise filtering}