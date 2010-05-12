\name{estep1}
\Rdversion{1.1}
\alias{estep1}
\title{
First E step of the AECM algorithm when fitting a mixture of PPCA models.
}
\description{
Internal function required for fitting a mixture of PPCA models.
}
\usage{
estep1(Y, Tau, Pi, mu, W, Sig, g, p, reset)
}
\arguments{
  \item{Y}{
A N x p data matrix.
}
  \item{Tau}{
A N x G matrix of posterior group membership probabilities.
}
  \item{Pi}{
A G vector of mixing proportions.
}
  \item{mu}{
A p x G matrix containing the mean for each group.
}
  \item{W}{
An p x q x G array of loadings for each group.
}
  \item{Sig}{
A scalar; the error covariance.
}
  \item{g}{
The number of groups currently being fitted.
}
  \item{p}{
Number of spectral bins in the NMR spectra.
}
  \item{reset}{
Logical indicating computational instability.
}
}
\details{
First E step of the AECM algorithm when fitting a mixture of PPCA models. An internal function.
}
\value{
A list containing
\item{Tau}{The N x G matrix of posterior group membership probablities.}
\item{logTau}{An N x G matrix of the log of the numerator of posterior group membership probablities.}
\item{reset}{Logical indicating computational instability.}
}
\references{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.
}
\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}
\note{
An internal function.
}

\seealso{
\code{\link{mppca.metabol}}
}

\keyword{internal}
