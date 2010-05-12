\name{mstep1}
\Rdversion{1.1}
\alias{mstep1}
\title{
First M-step of the AECM algorithm when fitting a mixture of PPCA models.
}
\description{
Internal function required for fitting a mixture of PPCA models.
}
\usage{
mstep1(Y, Tau, Pi, mu, g)
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
  \item{g}{
The number of groups currently being fitted.
}
}
\details{
First M-step of the AECM algorithm when fitting a mixture of PPCA models. An internal function.}
\value{
A list containing
\item{Pi}{A g vector of group probabilities}
\item{Mu}{A p x g matrix each column of which contains a group mean.}
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
