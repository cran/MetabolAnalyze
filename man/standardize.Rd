\name{standardize}
\Rdversion{1.1}
\alias{standardize}
\title{
Function to scale covariates.
}
\description{
A function to scale covariates so that they lie in [0,1] for reasons of stability and convergence of the EM algorithm.
}
\usage{
standardize(Covars)
}
\arguments{
  \item{Covars}{
An N x L matrix containing the L covariates of each of N observations.
}
}
\details{
A function to scale covariates so that they lie in [0,1] for reasons of stability and convergence of the EM algorithm. Care must be taken with categorical covariates: see \code{\link{ppcca.metabol}} for further information.
}
\value{
\item{Covars}{A standardized version of the input matrix of covariates.}
}

\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}


\keyword{internal }
