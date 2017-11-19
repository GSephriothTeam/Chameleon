\name{MixGOM}
\alias{MixGOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mixture Simulation based on generalized overlap of Maitra}
\description{
Generates a finite mixture model with Gaussian components for a prespecified level of goMega (generalized overlap of Maitra).
}
\usage{
MixGOM(goMega = NULL, K, p, sph = FALSE, hom = FALSE,
       ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), resN = 100,
       eps = 1e-06, lim = 1e06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{goMega}{ value of desired generalized overlap of Maitra. }
  \item{K}{ number of components. }
  \item{p}{ number of dimensions. }
  \item{sph}{ covariance matrix structure (FALSE - non-spherical, TRUE - spherical). }
  \item{hom}{ heterogeneous or homogeneous clusters (FALSE - heterogeneous, TRUE - homogeneous). }
  \item{ecc}{ maximum eccentricity. }
  \item{PiLow}{ value of the smallest mixing proportion (if 'PiLow' is not reachable with respect to K, equal proportions are taken; PiLow = 1.0 implies equal proportions by default). }
  \item{int}{ mean vectors are simulated uniformly on a hypercube with sides specified by int = (lower.bound, upper.bound). }
  \item{resN}{ maximum number of mixture resimulations. }
  \item{eps}{ error bound for overlap computation. }
  \item{lim}{ maximum number of integration terms (Davies, 1980). }
}
\details{
  Returns mixture parameters satisfying the prespecified level of goMega.
}
\value{
  \item{Pi }{vector of mixing proportions. }
  \item{Mu }{matrix consisting of components' mean vectors (K * p). }
  \item{S }{set of components' covariance matrices (p * p * K). }
  \item{goMega }{value of generalized overlap of Maitra. }
  \item{fail }{flag value; 0 represents successful mixture generation, 1 represents failure. }
}
\references{
Maitra, R. (2010) ``A re-defined and generalized percent-overlap-of-activation measure for studies of fMRI reproducibility and its use in identifying outlier activation maps'', NeuroImage, 50, 124-135.

Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.

Davies, R. (1980) ``The distribution of a linear combination of chi-square random variables'', Applied Statistics, 29, 323-333.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{ \code{overlapGOM}, \code{MixSim}, and \code{simdataset}. }
\examples{
%\dontrun{
set.seed(1234)

# controls average and maximum overlaps
(ex.1 <- MixGOM(goMega = 0.05, K = 4, p = 5))

# controls maximum overlap
(ex.2 <- MixGOM(goMega = 0.15, K = 4, p = 5, sph = TRUE))
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ datagen }% __ONLY ONE__ keyword per line
