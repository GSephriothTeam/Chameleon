\name{print.object}
\alias{print.MixSim}
\alias{summary.MixSim}
\alias{print.MixGOM}
\title{Functions for Printing or Summarizing Objects}
\description{
  A \code{MixSim} and \code{MixGOM} classes are declared, and these are functions
  to print and summarize objects.
}
\usage{
\method{print}{MixSim}(x, ...)
\method{summary}{MixSim}(object, ...)
\method{print}{MixGOM}(x, ...)
}
\arguments{
  \item{x}{an object with the 'MixSim' (or 'MixGOM') class attributes.}
  \item{object}{an object with the 'MixSim' (or 'MixGOM') class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  These are useful functions for summarizing and debugging.

  For other functions, they only show summaries of objects. Use
  \code{names} or \code{str} to explore the details.
}
\value{
  The results will cat or print on the STDOUT by default.
}
\references{
Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{
  \code{\link{MixSim}}.
}
\examples{
\dontrun{
# Functions applied by directly type the names of objects.
}
}
\keyword{programming}
