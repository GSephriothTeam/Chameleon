\name{overlap}
\alias{overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Overlap}
\description{
Computes misclassification probabilities and pairwise overlaps for finite mixture models with Gaussian components. Overlap is defined as sum of two misclassification probabilities.
}
\usage{
overlap(Pi, Mu, S, eps = 1e-06, lim = 1e06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pi}{ vector of mixing proprtions (length K). }
  \item{Mu}{ matrix consisting of components' mean vectors (K * p). }
  \item{S}{ set of components' covariance matrices (p * p * K). }
  \item{eps}{ error bound for overlap computation. }
  \item{lim}{ maximum number of integration terms (Davies, 1980). }
}

\value{
  \item{OmegaMap}{matrix of misclassification probabilities (K * K); OmegaMap[i,j] is the probability that X coming from the i-th component is classified to the j-th component.}
  \item{BarOmega}{value of average overlap.}
  \item{MaxOmega}{value of maximum overlap.}
  \item{rcMax}{row and column numbers for the pair of components producing maximum overlap 'MaxOmega'.}
}
\references{
Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.

Davies, R. (1980) ``The distribution of a linear combination of chi-square random variables'', Applied Statistics, 29, 323-333.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{ \code{MixSim}, \code{pdplot}, and \code{simdataset}. }
\examples{
%\dontrun{
data("iris", package = "datasets")
p <- ncol(iris) - 1
id <- as.integer(iris[, 5])
K <- max(id)

# estimate mixture parameters
Pi <- prop.table(tabulate(id))
Mu <- t(sapply(1:K, function(k){ colMeans(iris[id == k, -5]) }))
S <- sapply(1:K, function(k){ var(iris[id == k, -5]) })
dim(S) <- c(p, p, K)

overlap(Pi = Pi, Mu = Mu, S = S)
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

