\name{mstep2}
\Rdversion{1.1}
\alias{mstep2}
\title{
Second M-step of the AECM algorithm when fitting a mixture of PPCA models.
}
\description{
Internal function required for fitting a mixture of PPCA models.
}
\usage{
mstep2(Y, Tau, Pi, mu, W, Sig, g, p, q)
}
\arguments{
  \item{Y}{
A N x p data matrix.
}
  \item{Tau}{
An N x G matrix of posterior group membership probabilities.
}
  \item{Pi}{
A g vector of group probabilities.
}
  \item{mu}{
A p x g matrix containing the mean for each group.
}
  \item{W}{
A p x q x g array, each sheet of which contains a group specific loadings matrix.
}
  \item{Sig}{
The variance parmeter.
}
  \item{g}{
The number of groups currently being fitted.
}
  \item{p}{
The number of spectral bins in the NMR spectrum.
}
  \item{q}{
The number of principal components in the model being fitted.
}
}
\details{
Second M-step of the AECM algorithm when fitting a mixture of PPCA models. An internal function.
}
\value{
A list containing
\item{W}{A p x q x g array, each sheet of which contains a group specific loadings matrix.}
\item{Sig}{The variance parameter.}
}
\references{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.}
\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}
\note{
An internal function.
}

\seealso{
\code{\link{mppca.metabol}}
}
\keyword{ internal}
