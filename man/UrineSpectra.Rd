\name{UrineSpectra}
\Rdversion{1.1}
\alias{UrineSpectra}
\docType{data}
\title{
NMR metabolomic spectra from urine samples of 18 mice.
}
\description{
NMR metabolomic spectra from urine samples of 18 mice, each belonging to one of two treatment groups. Each spectrum has 189 spectral bins, measured in parts per million (ppm).

Covariates associated with the mice were also recorded: the weight of each mouse and the pH of their urine sample are provided. }
\usage{data(UrineSpectra)}
\format{
A list containing
\enumerate{
\item a matrix with 18 rows and 189 columns
\item a data frame with 18 observations on 3 variables:
\itemize{
  \item Treatment group membership of each animal.
  \item Weight (in grammes) of each animal.       
  \item pH of the urine sample from each animal.
}
}
}
\details{
This is simulated data, based on parameter estimates from a PPCA model with two prinicipal components fitted to a similar real data set.
}
\source{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.
}


\keyword{datasets}
