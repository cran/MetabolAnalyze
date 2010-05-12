\name{BrainSpectra}
\Rdversion{1.1}
\alias{BrainSpectra}
\docType{data}
\title{
NMR spectral data from brain tissue samples.
}
\description{
NMR spectral data from brain tissue samples of 33 rats, where each tissue sample originates in one of four known brain regions. Each spectrum has 164 spectral bins, measured in parts per million (ppm).

}
\usage{data(BrainSpectra)}
\format{
  A list containing
\enumerate{
\item a matrix with 33 rows and 164 columns
\item a vector indicating the brain region of origin of each sample where:
\itemize{
  \item 1 = Brain stem
  \item 2 = Cerebellum
  \item 3 = Hippocampus
  \item 4 = Pre-frontal cortex
}
}}
\details{
This is simulated data, based on parameter estimates from a mixture of PPCA models with 4 groups and 7 principal components fitted to a similar real data set.
}
\source{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.
}

\keyword{datasets}
