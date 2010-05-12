\name{mppca.loadings.plot}
\Rdversion{1.1}
\alias{mppca.loadings.plot}
\title{
Plot loadings resulting from fitting a MPPCA model.
}
\description{
A function to plot the loadings resulting from fitting a MPPCA model to metabolomic data. A barplot or a scatterplot can be produced.
}
\usage{
mppca.loadings.plot(output, Y, barplot = FALSE, labelsize = 0.3)
}
\arguments{
  \item{output}{
An object resulting from fitting a MPPCA model.
}
  \item{Y}{
The N x p matrix of observations to which the MPPCA model is fitted.
}
  \item{barplot}{
Logical indicating whether a barplot of the loadings is required rather than a scatter plot. By default a scatter plot is produced.
}
  \item{labelsize}{
Size of the text of the spectral bin labels on the resulting plot.
}
}
\details{
A function which produces a series of plots illustrating the loadings resulting from fitting a MPPCA model to metabolomic data. A barplot or a scatterplot can be produced. The size of the text of the spectral bin labels on the bar plot can also be adjusted if the number of bins plotted is large.
}

\references{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.}

\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}

\seealso{
\code{\link{mppca.metabol}}
}


\keyword{ methods }
