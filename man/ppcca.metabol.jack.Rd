\name{ppcca.metabol.jack}
\Rdversion{1.1}
\alias{ppcca.metabol.jack}
\title{
Fit a probabilistic principal components and covariates analysis model to a metabolomic data set, and assess uncertainty via the jackknife.
}
\description{
Fit a probabilistic principal components and covariates analysis (PPCCA) model to a metabolomic data set via the EM algorithm, and assess uncertainty in the obtained loadings estimates and the regression coefficients via the jackknife.
}
\usage{
ppcca.metabol.jack(Y, Covars, minq=1, maxq=2, scale="none", epsilon=0.1, 
conflevel=0.95)
}

\arguments{
  \item{Y}{
An N x p data matrix in which each row is a spectrum.
}
  \item{Covars}{
An N x L covariate data matrix where each row is a set of covariates.
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
Level of confidence required for the loadings and regression coefficients confidence intervals. By default 95\eqn{\%} confidence intervals are computed.
}
}
\details{
A (range of) PPCCA model(s) are fitted and an optimal model (i.e. number of principal components, q) is selected. Confidence intervals for the obtained loadings and regression coefficients are then obtained via the jackknife i.e. a model with q principal components is fitted to the data \eqn{N} times, where an observation is removed from the dataset each time. 

Care should be taken with the form of covariates supplied. All covariates are standardized (to lie in [0,1]) within the ppcca.metabol.jack function for stability reasons. Hence continuous covariates and binary valued categorical covariates are easily handled. For a categorical covariate with V levels, the equivalent V-1 dummy variables representation should be passed as an argument to ppcca.metabol.jack.

On convergence of the algorithm, the number of loadings significantly different from zero is printed on screen. The user may then further examine the significant loadings when prompted by selecting a cutoff value from the table printed on screen. Bar plots detailing the resulting significantly high loadings are provided.
}
\value{
A list containing:
\item{q}{The number of principal components in the optimal PPCCA model, selected by the BIC.}
\item{sig}{The posterior mode estimate of the variance of the error terms.}
\item{scores}{An N x q matrix of estimates of the latent locations of each observation in the principal subspace.}
\item{loadings}{The maximum likelihood estimate of the p x q loadings matrix.}
\item{SignifW}{The maximum likelihood estimate of the loadings matrix for those loadings significantly different from zero.}
\item{SignifHighW}{The maximum likelihood estimate of the loadings matrix for those loadings significantly different from zero and above the user selected cutoff point.}
\item{LowerCI_W}{The lower limit of the confidence interval for those loadings significantly different from zero.}
\item{UpperCI_W}{The upper limit of the confidence interval for those loadings significantly different from zero.}
\item{coefficients}{The maximum likelihood estimates of the regression coefficients.}
\item{coeffCI}{A matrix detailing the upper and lower limits of the confidence intervals for the regression parameters.}
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
\code{\link{ppcca.metabol}}, \code{\link{ppcca.scores.plot}},\code{\link{loadings.jack.plot}} 
}
\examples{
data(UrineSpectra)
\dontrun{
mdlfit<-ppcca.metabol.jack(UrineSpectra[[1]], UrineSpectra[[2]][,2], minq=2, maxq=2)
loadings.jack.plot(mdlfit)
ppcca.scores.plot(mdlfit, UrineSpectra[[2]][,2], group=UrineSpectra[[2]][,1], covarnames="Weight")
}
}
\keyword{ methods }

