\name{GOdata}

\alias{GOdata}
\alias{resultFisher}
\alias{resultKS}

\docType{data}

\title{Sample topGOdata and topGOresult objects}

\description{
  The \code{GOdata} contains an instance of a \code{topGOdata} object.
  It can be used to run an enrichment analysis directly.

  The \code{resultFisher} contains the results of an enrichment
  analysis. 
}

\usage{data(GOdata)}

\source{
  Generated using the \code{ALL} gene expression data. See
  \code{\link{topGOdata-class}} for code examples on how-to generate
  such an object.
}

\examples{
data(GOdata)

## print the object
GOdata

data(results.tGO)

## print the object
resultFisher
}

\keyword{datasets}