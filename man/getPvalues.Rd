\name{getPvalues}

\alias{getPvalues}


\title{Convenient function to compute p-values from a gene expression matrix.}
\description{Warping function of "mt.teststat", for computing p-values
  of a gene expression matrix.
}

\usage{
   getPvalues(edata, classlabel, test = "t", 
	      alternative = c("greater", "two.sided", "less")[1], genesID = NULL, 
	      correction = c("none", "Bonferroni", "Holm", "Hochberg", "SidakSS", 
			     "SidakSD", "BH", "BY")[8]) 
}

\arguments{
  \item{edata}{Gene expression matrix.}
  \item{classlabel}{The phenotype of the data}
  \item{test}{Which test statistic to use}
  \item{alternative}{The alternative of the test statistic}
  \item{genesID}{if a subset of genes is provided}
  \item{correction}{Multiple testing correction procedure}
}

\value{
  An named numeric vector of p-values.
}

\author{Adrian Alexa}


\seealso{
  \code{\link{GOKSTest}},  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}


\examples{

library(ALL)
data(ALL)

## discriminate B-cell from T-cell
classLabel <- as.integer(sapply(ALL$BT, function(x) return(substr(x, 1, 1) == 'T')))

## Differentially expressed genes
geneList <- getPvalues(exprs(ALL), classlabel = classLabel,
                       alternative = "greater", correction = "BY")

hist(geneList, 50)
}

\keyword{graphs}
