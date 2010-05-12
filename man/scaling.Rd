\name{scaling}
\Rdversion{1.1}
\alias{scaling}
\title{
Function to scale metabolomic spectral data.
}
\description{
This function provides the options of Pareto scaling, unit scaling or no scaling of metabolomic data.
}
\usage{
scaling(Y, type = "none")
}
\arguments{
  \item{Y}{
  	An N x p matrix of metabolomic spectra. Each row of Y is an observation's spectrum.
}
  \item{type}{
Default is "none" meaning the data are not altered. If "pareto", the data are Pareto scaled. If "unit", the data are unit scaled.
}
}
\details{
Pareto scaling, frequently utilised in metabolomic analyses, scales data by dividing each variable by the square root of the standard deviation. Unit scaling divides each variable by the standard deviation so that each variable has variance equal to 1.
}
\value{
The function returns the requested scaled version of the input matrix \eqn{Y}.
}
\references{
van den Berg, R.A., Hoefsloot, H.C.J, Westerhuis, J.A. and Smilde, A.K. and van der Werf, M.J. (2006) Centering, scaling, and transformations: improving the biological information content of metabolomics data. \emph{BMC Genomics}, 7, 1, 142.

Nyamundanda G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report, University College Dublin.
}
\author{
Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan.
}


\keyword{ internal }
