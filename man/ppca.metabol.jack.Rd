\name{ppca.metabol.jack}
\Rdversion{1.1}
\alias{ppca.metabol.jack}
\title{
Fit a probabilistic principal components analysis model to a metabolomic data set, and assess uncertainty via the jackknife.
}
\description{
Fit a probabilistic principal components analysis (PPCA) model to a metabolomic data set via the EM algorithm, and assess uncertainty in the obtained loadings estimates via the jackknife.
}
\usage{
ppca.metabol.jack(Y, minq=1, maxq=2, scale ="none", 
epsilon = 0.1, conflevel = 0.95)
}
\arguments{
  \item{Y}{
An N x p data matrix where each row is a spectrum.
}
  \item{minq}{
The minimum number of principal components to be fit. By default minq is 1.
}
  \item{maxq}{
The maximum number of principal components to be fit. By default maxq is 2.
}
  \item{scale}{
Type of scaling of the data which is required. The default is "none". Options include "pareto' and "unit" scaling. See \code{\link{scaling}} for further details.
}
  \item{epsilon}{
Value on which the convergence assessment criterion is based. Set by default to 0.1.
}
  \item{conflevel}{
Level of confidence required for the loadings confidence intervals. By default 95\eqn{\%} confidence intervals are computed.
}
}
\details{
A (range of) PPCA model(s) are fitted and an optimal model (i.e. number of principal components, q) is selected. Confidence intervals for the obtained loadings are then obtained via the jackknife i.e. a model with q principal components is fitted to the dataset \eqn{N} times, where an observation is removed from the dataset each time. 

On convergence of the algorithm, the number of loadings significantly different from zero is printed on screen. The user may then further examine the significant loadings when prompted by selecting a cutoff value from the table printed on screen. Bar plots detailing the resulting significantly high loadings are provided.
}
\value{
A list containing:
\item{q}{The number of principal components in the optimal PPCA model, selected by the BIC.}
\item{sig}{The maximum likelihood estimate of the variance of the error terms.}
\item{scores}{An N x q matrix of estimates of the latent locations of each observation in the principal subspace.}
\item{loadings}{The maximum likelihood estimate of the p x q loadings matrix.}
\item{SignifW}{The maximum likelihood estimate of the loadings matrix for those loadings significantly different from zero.}
\item{SignifHighW}{The maximum likelihood estimate of the loadings matrix for those loadings significantly different from zero and higher than a user selected cutoff point.}
\item{Lower}{The lower limit of the confidence interval for those loadings significantly different from zero.}
\item{Upper}{The upper limit of the confidence interval for those loadings significantly different from zero.}
\item{Cutoffs}{A table detailing a range of cutoff points and the associated number of selected spectral bins.}
\item{number}{The number of spectral bins selected by the user.}
\item{cutoff}{The cutoff value selected by the user.}
\item{BIC}{A vector containing the BIC values for the fitted models.}
\item{AIC}{A vector containing the AIC values for the fitted models.}
}
\references{
Nyamundanda G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report, University College Dublin.

}
\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan.
}

\seealso{
\code{\link{ppca.metabol}}, \code{\link{loadings.jack.plot}}, \code{\link{ppca.scores.plot}}
}
\examples{
data(UrineSpectra)
\dontrun{
mdlfit<-ppca.metabol.jack(UrineSpectra[[1]], minq=2, maxq=2, scale="none")
loadings.jack.plot(mdlfit)
ppca.scores.plot(mdlfit, group=UrineSpectra[[2]][,1])}
}

\keyword{ methods }
