\name{MixSim-package}
\alias{00_MixSim-package}
\docType{package}
\title{
Simulation of Gaussian Finite Mixture Models
}
\description{
Simulation of Gaussian finite mixture models for prespecified levels of average or/and maximum overlap. Pairwise overlap is defined as the sum of two misclassification probabilities.
}
\details{
\tabular{ll}{
Package: \tab MixSim\cr
Type: \tab Package\cr
% Version: \tab 1.0\cr
Date: \tab 2012-08-12\cr
License: \tab GPL (>= 2)\cr
}
Function 'MixSim' simulates a finite mixture model for a prespecified level of average or/and maximum overlap.

Function 'overlap' computes all misclassification probabilities for a finite mixture model.

Function 'pdplot' constructs a parallel distribution plot for a finite mixture model.

Function 'simdataset' simulates a dataset from a finite mixture model.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra.

Maintainer: Volodymyr Melnykov <vmelnykov@cba.ua.edu>
}
\references{
Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.

Davies, R. (1980) ``The distribution of a linear combination of chi-square random variables'', Applied Statistics, 29, 323-333.

Meila, M. (2006) ``Comparing clusterings - an information based distance'', Journal of Multivariate Analysis, 98, 873-895.
}
\examples{
%\dontrun{
# Simulate parameters of a mixture model
A <- MixSim(BarOmega = 0.01, MaxOmega = 0.10, K = 10, p = 5)

# Display the mixture via the parallel distribution plot
pdplot(A$Pi, A$Mu, A$S, MaxInt = 0.5)
%}
}
