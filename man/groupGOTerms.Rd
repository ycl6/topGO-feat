\name{groupGOTerms}
\alias{groupGOTerms}

\alias{GOCCTerm}
\alias{GOBPTerm}
\alias{GOMFTerm}


\title{Grouping of GO terms into the three ontologies}
\description{
  This function split the GOTERM environment into three different
  ontologies. The newly created environments contain each only
  the terms from one of the following ontologies 'BP', 'CC', 'MF'
}

\usage{
groupGOTerms(where)
}

\arguments{
  \item{where}{The the environment where you want to bind the results.}
}

\value{
  The function returns NULL.
}


\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link[GO.db:GOTERM]{GOTERM}}
}

\examples{
groupGOTerms()
}

\keyword{misc}
