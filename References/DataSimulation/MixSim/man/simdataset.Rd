\name{simdataset}
\alias{simdataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dataset Simulation }
\description{
  Simulates a datasets of sample size n given parameters of finite mixture model with Gaussian components.
}
\usage{
simdataset(n, Pi, Mu, S, n.noise = 0, n.out = 0, alpha = 0.001,
           max.out = 100000, int = NULL, lambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ sample size. }
  \item{Pi}{ vector of mixing proportions (length K). }
  \item{Mu}{ matrix consisting of components' mean vectors (K * p). }
  \item{S}{ set of components' covariance matrices (p * p * K). }
  \item{n.noise}{ number of noise variables. }
  \item{n.out}{ number of outlying observations. }
  \item{alpha}{ level for simulating outliers. }
  \item{max.out}{ maximum number of trials to simulate outliers. }
  \item{int}{ interval for noise and outlier generation. }
  \item{lambda}{ inverse Box-Cox transformation coefficients. }
}
\details{
   The function simulates a dataset of n observations from a mixture model with parameters 'Pi' (mixing proportions), 'Mu' (mean vectors), and 'S' (covariance matrices). Mixture component sample sizes are produced as a realization from a multinomial distribution with probabilities given by mixing proportions.
   To make a dataset more challenging for clustering, a user might want to simulate noise variables or outliers. Parameter 'n.noise' specifies the desired number of noise variables. If an interval 'int' is specified, noise will be simulated from a Uniform distribution on the interval given by 'int'. Otherwise, noise will be simulated uniformly between the smallest and largest coordinates of mean vectors. 'n.out' specifies the number of observations outside (1 - 'alpha') ellipsoidal contours for the weighted component distributions. Outliers are simulated on a hypercube specified by the interval 'int'.
   A user can apply an inverse Box-Cox transformation providing a vector of coefficients 'lambda'. The value 1 implies that no transformation is needed for the corresponding coordinate.
}
\value{
  \item{X }{ simulated dataset (n + n.out) x (p + n.noise); noise coordinates are provided in the last n.noise columns.}
  \item{id }{ classification vector (length n + n.out); 0 represents an outlier.}
}
\references{
Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{ \code{MixSim}, \code{overlap}, and \code{pdplot}. }
\examples{
\dontrun{
set.seed(1234)

repeat{
   Q <- MixSim(BarOmega = 0.01, K = 4, p = 2)
   if (Q$fail == 0) break
}

# simulate a dataset of size 300 and add 10 outliers simulated on (0,1)x(0,1)
A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.out = 10, int = c(0, 1))
colors <- c("red", "green", "blue", "brown", "magenta")
plot(A$X, xlab = "x1", ylab = "x2", type = "n")
for (k in 0:4){
   points(A$X[A$id == k, ], col = colors[k+1], pch = 19, cex = 0.5)
}

repeat{
   Q <- MixSim(MaxOmega = 0.1, K = 4, p = 1)
   if (Q$fail == 0) break
}

# simulate a dataset of size 300 with 1 noise variable
A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.noise = 1)
plot(A$X, xlab = "x1", ylab = "x2", type = "n")
for (k in 1:4){
   points(A$X[A$id == k, ], col = colors[k+1], pch = 19, cex = 0.5)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ datagen }% __ONLY ONE__ keyword per line
