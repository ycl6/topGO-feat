\name{geneList}

\alias{geneList}
\alias{topDiffGenes}
\alias{affyLib}


\docType{data}

\title{A toy example of a list of gene identifiers and the respective p-values}
\description{
  The \code{geneList} data is compiled from a differential expression
  analysis of the \code{ALL} dataset. It contains just a small number
  of genes with the corespondent p-values. The information on where
  to find the GO annotations is stored in the \code{ALL} object.
  
  The \code{topDiffGenes} function included in this dataset will select 
  the differentially expressed genes, at 0.01 significance level, from
  \code{geneList}.
}

\usage{data(geneList)}

\source{
  Generated using the \code{ALL} gene expression data. See the "scripts" directory. 
}

\examples{
data(geneList)

## print the object
head(geneList)
length(geneList)

## the number of genes with a p-value less than 0.01
sum(topDiffGenes(geneList))
}

\keyword{datasets}