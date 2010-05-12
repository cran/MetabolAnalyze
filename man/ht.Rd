\name{ht}
\Rdversion{1.1}
\alias{ht}
\title{
Function to plot a heatmap of BIC values.
}
\description{
Function to plot a heat map of BIC values where lighter colours indicate larger values and optimal models. A black cross indicates the optimal model.

The function is a modified version of \code{\link{heatmap}}.
}
\usage{
ht(x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, distfun = dist, hclustfun = hclust, reorderfun = function(d, w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, "Rowv"), scale = c("row", "column", "none"), na.rm = FALSE, margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 1, cexCol = 1, labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, verbose = getOption("verbose"), q, g)
}
\arguments{
See the help file for \code{\link{heatmap}}.
}
\details{
This function is used internally in \code{\link{mppca.metabol}}.
}
\value{
See the help file for \code{\link{heatmap}}.
}
\references{
Nyamundanda, G., Gormley, I.C. and Brennan, L. (2010) Probabilistic principal components analysis for metabolomic data. Technical report. University College Dublin, Ireland.
}
\author{
	Nyamundanda Gift, Isobel Claire Gormley and Lorraine Brennan
}
\note{
An internal function.
}

\seealso{
\code{\link{mppca.metabol}}
}

\keyword{ internal }

