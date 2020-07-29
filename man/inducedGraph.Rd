\name{inducedGraph}

\alias{inducedGraph}
\alias{nodesInInducedGraph}


\title{The subgraph induced by a set of nodes.}
\description{
  Given a set of nodes (GO terms) this function is returning
  the subgraph containing these nodes and their ancestors.
}

\usage{
inducedGraph(dag, startNodes)
nodesInInducedGraph(dag, startNodes)
}

\arguments{
  \item{dag}{An object of class \code{graphNEL} containing a directed
    graph.}
  \item{startNodes}{A character vector giving the starting nodes.}
}

\value{
  An object of class \code{\link[graph]{graphNEL-class}} is returned.
}

\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link{reverseArch}},
}

\examples{
data(GOdata)

## the GO graph
g <- graph(GOdata)
g

## select 10 random nodes
sn <- sample(nodes(g), 10)


## the subgraph induced by these nodes
sg <- inducedGraph(g, sn)
sg
}

\keyword{graphs}
