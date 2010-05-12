\name{ppcca.scores.plot}
\Rdversion{1.1}
\alias{ppcca.scores.plot}
\title{
Plot scores from a fitted PPCCA model.
}
\description{
A function to plot the scores resulting from fitting a PPCCA model to metabolomic data.
}
\usage{
ppcca.scores.plot(output, Covars, group = FALSE)
}
\arguments{
  \item{output}{
An object resulting from fitting a PPCCA model.
}
  \item{Covars}{
An N x L covariate data matrix where each row is a set of covariates.
}
  \item{group}{
Should it be relevant, a vector indicating the known treatment group membership of each observation.
}
}
\details{
This function produces a series of scatterplots each illustrating the estimated score for each observation within the reduced q dimensional space. The uncertainty associated with the score estimate is also illustrated through its 95% posterior set. Additionally, the size of the plotted point for each observation is proportional to the size of the covariate for that observation.

It is often the case that observations are known to belong to treatment groups; the treatment group membership of each observation can be illustrated on the plots produced by utilizing the `group' argument.
}

\references{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.}

\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}


\seealso{
\code{\link{ppcca.metabol}}, \code{\link{ppcca.metabol.jack}}
}



\keyword{ methods }