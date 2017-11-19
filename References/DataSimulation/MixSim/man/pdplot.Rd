\name{pdplot}
\alias{pdplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallel Distribution Plot }
\description{
Constructs a parallel distribution plot for Gaussian finite mixture models.
}
\usage{
pdplot(Pi, Mu, S, file = NULL, Nx = 5, Ny = 5, MaxInt = 1, marg = c(2,1,1,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pi}{ vector of mixing proportions. }
  \item{Mu}{ matrix consisting of components' mean vectors (K * p). }
  \item{S}{ set of components' covariance matrices (p * p * K). }
  \item{file}{ name of .pdf-file. }
  \item{Nx}{ number of color levels for smoothing along the x-axis. }
  \item{Ny}{ number of color levels for smoothing along the y-axis. }
  \item{MaxInt}{ maximum color intensity. }
  \item{marg}{ plot margins. }
}
\details{
If 'file' is specified, produced plot will be saved as a .pdf-file.
}

\references{
Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{ \code{MixSim}, \code{overlap}, and \code{simdataset}. }
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

pdplot(Pi = Pi, Mu = Mu, S = S)
%}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ hplot }% __ONLY ONE__ keyword per line
