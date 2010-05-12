\name{loadings.jack.plot}
\Rdversion{1.1}
\alias{loadings.jack.plot}
\title{
Plot loadings and their associated confidence intervals.
}
\description{
A function to plot the loadings and confidence intervals resulting from fitting a PPCA model or a PPCCA model to metabolomic data.
}
\usage{
loadings.jack.plot(output)
}
\arguments{
  \item{output}{
An object resulting from fitting a PPCA model or a PPCCA model.
}
}
\details{
The function produces a plot of those loadings on the first principal component which are significantly different from zero, and higher than a user specified cutoff point. Error bars associated with the estimates, derived using the jackknife, are also plotted.
}

\references{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.}

\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}

\seealso{
\code{\link{ppca.metabol.jack}}, \code{\link{ppcca.metabol.jack}}
}


\keyword{methods}
