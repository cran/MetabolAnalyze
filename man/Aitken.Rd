\name{Aitken}
\Rdversion{1.1}
\alias{Aitken}
\title{
Assess convergence of an EM algorithm.
}
\description{
This function assesses convergence of the EM algorithm using Aitken's acceleration method, when fitting a PPCA based model.
}
\usage{
Aitken(ll, lla, v, q, epsilon)
}
\arguments{
  \item{ll}{
A vector of log likelihoods from the current and previous iterations.
}
  \item{lla}{
A vector containing the asympototic estimates of the maximized log likelihoods from the current and previous iterations.
}
  \item{v}{
Iteration number.
}
  \item{q}{
The dimension of the latent principal subspace for the PPCA based model currently being fitted.
}
\item{epsilon}{
The value on which convergence of the EM algorithm is based. 
}
}
\details{
This function assesses convergence of the EM algorithm using Aitken's acceleration method in which an estimate of the maximized log likelihood at each iteration is evaluated. Convergence is achieved when the absolute difference between contiguous estimates, tol, is less than some user defined level, epsilon.
}
\value{
A list containing:
\item{tol }{The absolute difference between contiguous estimates of the asymptotic maximized log likelihood.}
\item{la }{The asymptotic estimate of the maximized log likelihood at the current iteration.}
}
\references{
McLachlan, G.J. and Krishnan, T. (1997) The EM algorithm and Extensions. Wiley, New York.}
\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}
\note{
This is used internally in functions which fit PPCA based models via the EM algorithm within the package MetabolAnalyze.
}
\seealso{
\code{\link{ppca.metabol}}, \code{\link{ppcca.metabol}}, \code{\link{mppca.metabol}}
}

\keyword{ internal }


