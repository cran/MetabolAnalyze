\name{ppcca.metabol}
\Rdversion{1.1}
\alias{ppcca.metabol}
\title{
Fit a probabilistic principal components and covariates analysis (PPCCA) model to a metabolomic data set via the EM algorithm.
}
\description{
This function fits a probabilistic principal components and covariates analysis model to metabolomic spectral data via the EM algorithm.
}
\usage{
ppcca.metabol(Y, Covars, minq=1, maxq=2, scale = "none", epsilon = 0.1, 
plot.BIC = FALSE, printout=TRUE)
}
\arguments{
  \item{Y}{
An N x p data matrix in which each row is a spectrum.
}
  \item{Covars}{
An N x L covariate data matrix in which each row is a set of covariates.
}
  \item{minq}{
The minimum number of principal components to be fit.
}
  \item{maxq}{
The maximum number of principal components to be fit.
}
  \item{scale}{
Type of scaling of the data which is required. The default is "none". Options include "pareto' and "unit" scaling. See \code{\link{scaling}} for further details.
}
  \item{epsilon}{
Value on which the convergence assessment criterion is based. Set by default to 0.1.
}
  \item{plot.BIC}{
Logical indicating whether or not a plot of the BIC values for the different models fitted should be provided. By default, the plot is not produced.
}
\item{printout}{
Logical indicating whether or not a statement is printed on screen detailing the progress of the algorithm.
}
}
\details{
This function fits a probabilistic principal components and covariates analysis model to metabolomic spectral data via the EM algorithm. A range of models with different numbers of principal components can be fitted.

Care should be taken with the form of covariates supplied. All covariates are standardized (to lie in [0,1]) within the ppcca.metabol function for stability reasons. Hence continuous covariates and binary valued categorical covariates are easily handled. For a categorical covariate with V levels, the equivalent V-1 dummy variables representation should be passed as an argument to ppcca.metabol.
}
\value{
A list containing:
\item{q}{The number of principal components in the optimal PPCCA model, selected by the BIC.}
\item{sig}{The maximum likelihood estimate of the variance of the error terms.}
\item{scores}{An N x q matrix of estimates of the latent locations of each observation in the principal subspace.}
\item{loadings}{The maximum likelihood estimate of the p x q loadings matrix.}
\item{coefficients}{The maximum likelihood estimates of the regression coefficients associated with the covariates in the PPCCA model.}
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
\code{\link{ppcca.metabol.jack}}, \code{\link{ppcca.scores.plot}}\code{\link{loadings.plot}}
}
\examples{
data(UrineSpectra)
mdlfit<-ppcca.metabol(UrineSpectra[[1]], UrineSpectra[[2]][,2:3], minq=2, maxq=2)
loadings.plot(mdlfit)
ppcca.scores.plot(mdlfit, UrineSpectra[[2]][,2:3], UrineSpectra[[2]][,1])
}


\keyword{ methods }

