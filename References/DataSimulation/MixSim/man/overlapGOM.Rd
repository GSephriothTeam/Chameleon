\name{overlapGOM}
\alias{overlapGOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized overlap of Maitra}
\description{
Computes the generalized overlap as defined by R. Maitra.
}
\usage{
overlapGOM(Pi, Mu, S, eps = 1e-06, lim = 1e06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pi}{ vector of mixing proprtions (length K). }
  \item{Mu}{ matrix consisting of components' mean vectors (K * p). }
  \item{S}{ set of components' covariance matrices (p * p * K). }
  \item{eps}{ error bound for overlap computation. }
  \item{lim}{ maximum number of integration terms (Davies, 1980). }
}

\value{Returns the value of goMega.}
\references{
Maitra, R. (2010) ``A re-defined and generalized percent-overlap-of-activation measure for studies of fMRI reproducibility and its use in identifying outlier activation maps'', NeuroImage, 50, 124-135.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.

Davies, R. (1980) ``The distribution of a linear combination of chi-square random variables'', Applied Statistics, 29, 323-333.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{ \code{MixSim}, \code{MixGOM}, and \code{overlap}. }
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

overlapGOM(Pi = Pi, Mu = Mu, S = S)
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

