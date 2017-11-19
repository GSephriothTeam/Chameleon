\name{VarInf}
\alias{VarInf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variation of Information}
\description{
Computes the variation of information for two classification vectors.
}
\usage{
VarInf(id1, id2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id1}{ first partitioning vector. }
  \item{id2}{ second partitioning vector. }
}

\value{
  {Returns the variation of information. It is equal to 0 if and only if two classification vectors are identical.}
}
\references{
Meila, M. (2006) ``Comparing clusterings - an information based distance'', Journal of Multivariate Analysis, 98, 873-895.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{ \code{ClassProp}, and \code{RandIndex}. }
\examples{
%\dontrun{
id1 <- c(rep(1, 50), rep(2,100))
id2 <- rep(1:3, each = 50)
VarInf(id1, id2)
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

