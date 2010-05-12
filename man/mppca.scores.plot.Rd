\name{mppca.scores.plot}
\Rdversion{1.1}
\alias{mppca.scores.plot}
\title{
Plot scores from a fitted MPPCA model
}
\description{
A function to plot the scores resulting from fitting a MPPCA model to metabolomic data.
}
\usage{
mppca.scores.plot(output, group = FALSE, gplegend = TRUE)
}
\arguments{
  \item{output}{
An object resulting from fitting a MPPCA model.
}
  \item{group}{
Should it be relevant, a vector indicating the known treatment group membership of each observation prior to clustering.
}
  \item{gplegend}{
Logical indicating whether a legend should be plotted.
}
}
\details{
This function produces a series of scatterplots, for each group uncovered. For group g, each scatterplot illustrates the estimated score for each observation allocated to that group within the reduced q dimensional space. The uncertainty associated with the score estimate is also illustrated through its 95% posterior set. 

It is often the case that observations are known to belong to treatment groups, for example, and the MPPCA model is employed to uncover any underlying subgroups, possibly related to disease subtypes. The treatment group membership of each observation can be illustrated on the plots produced by utilizing the `group' argument.
}

\references{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.}

\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}


\seealso{
\code{\link{mppca.metabol}}
}



\keyword{ methods}
