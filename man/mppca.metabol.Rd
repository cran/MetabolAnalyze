\name{mppca.metabol}
\Rdversion{1.1}
\alias{mppca.metabol}
\title{
Fit a mixture of probabilistic principal components analysis (MPPCA) model to a metabolomic data set via the EM algorithm to perform simultaneous dimension reduction and clustering.
}
\description{
This function fits a mixture of probabilistic principal components analysis model to metabolomic spectral data via the EM algorithm.
}
\usage{
mppca.metabol(Y, minq=1, maxq=2, ming, maxg, scale = "none", 
epsilon = 0.1, plot.BIC = FALSE)
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
  \item{ming}{
The minimum number of groups to be fit.
}
  \item{maxg}{
The maximum number of groups to be fit.
}
  \item{scale}{
Type of scaling of the data which is required. The default is "none". Options include "pareto' and "unit" scaling. See \code{\link{scaling}} for further details.
}
  \item{epsilon}{
Value on which the convergence assessment criterion is based. Set by default to 0.1.
}
  \item{plot.BIC}{
Logical indicating whether or not a plot of the BIC values for the different models fitted should be provided. By default, the plot is not produced.}
}
\details{
This function fits a mixture of probabilistic principal components analysis models to metabolomic spectral data via the EM algorithm. A range of models with different numbers of groups and different numbers of principal components can be fitted. The model performs simultaneous clustering of observations into unknown groups and dimension reduction simultaneously.
}
\value{
A list containing:
\item{q}{The number of principal components in the optimal MPPCA model, selected by the BIC.}
\item{g}{The number of groups in the optimal MPPCA model, selected by the BIC.}
\item{sig}{The posterior mode estimate of the variance of the error terms.}
\item{scores}{A list of length g, each entry of which is a n_g x q matrix of estimates of the latent locations of each observation in group g in the principal subspace.}
\item{loadings}{An array of dimension p x q x g, each sheet of which contains the maximum likelihood estimate of the p x q loadings matrix for a group.}
\item{Pi}{The vector indicating the probability of belonging to each group.}
\item{mean}{A p x g matrix, each column of which contains a group mean.}
\item{tau}{An N x g matrix, each row of which contains the posterior group membership probabilities for an observation.}
\item{clustering}{A vector of length N indicating the group to which each observation belongs.}
\item{BIC}{A matrix containing the BIC values for the fitted models.}
\item{AIC}{A matrix containing the AIC values for the fitted models.}
}
\references{
Nyamundanda G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report, University College Dublin.
}
\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan.
}

\seealso{
\code{\link{mppca.scores.plot}}, \code{\link{mppca.loadings.plot}}
}
\examples{
data(BrainSpectra)
\dontrun{
mdlfit<-mppca.metabol(BrainSpectra[[1]], minq=7, maxq=7, ming=4, maxg=4, 
plot.BIC = TRUE)
mppca.scores.plot(mdlfit)
mppca.loadings.plot(mdlfit, BrainSpectra[[1]])
}
}

\keyword{methods}
