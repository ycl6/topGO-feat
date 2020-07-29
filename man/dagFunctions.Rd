\name{Determines the levels of a Directed Acyclic Graph (DAG)}

\alias{buildLevels}
\alias{getNoOfLevels}
\alias{getGraphRoot}
\alias{reverseArch}



\title{Utility functions to work with Directed Acyclic Graphs (DAG)}

\description{Basic functions to work witg DAGs}

\usage{
buildLevels(dag, root = NULL, leafs2root = TRUE)
getNoOfLevels(graphLevels)
getGraphRoot(dag, leafs2root = TRUE)
reverseArch(dirGraph, useAlgo = "sparse", useWeights = TRUE)
}


\arguments{
  \item{dag}{A \code{graphNEL} object.}
  
  \item{root}{A character vector specifing the root(s) of the DAG. If
    not specified the root node is autmatically computed.}

  \item{leafs2root}{The leafs2root parameter tell if the graph has edges
    directed from the leaves to the root, or vice-versa}

  \item{graphLevels}{An object of type list, returned by the \code{buildLevels} function.}  
  
  \item{dirGraph}{A \code{graphNEL} object containing a directed graph.}
  
  \item{useAlgo}{A character string specifing one of the following options
    \code{c("sparse", "normal")}. By default, \code{useAlgo = "sparse"}, a sparce matrix object is
    used to transpose the adjacency matrix. Otherwise a standard R
    martix is used.}
  
  \item{useWeights}{If weights should be used (if \code{useAlgo =
      "normal"} then the weigths are used anyway)}
}


\details{
  \code{buildLevels} function determines the levels of a Directed
  Acyclic Graph (DAG). The level of a node is defined as the longest
  path from the node to the root. The function take constructs a named
  list containg varios information about each nodes level. The root has
  level 1.

  \code{getNoOfLevels} - a convenient function to extract the number of
  levels from the object returned by \code{buildLevels}
  
  \code{getGraphRoot} finds the root(s) of the DAG

  \code{reverseArch} - simple function to invert the direction of edges
  in a DAG. The returned graph is of class graphNEL.  It can use either
  simple matrices or sparse matrices (SparseM library)
}

\value{
  \code{buildLevels} returns a list containing:
  \item{level2nodes}{Environment where the key is the level number
    with the value being the nodes on that level.}
  \item{nodes2level}{Environment where the key is the node label (the GO
    ID) and the value is the level on which that node lies.}
  \item{noOfLevels}{The number of levels}
  \item{noOfNodes}{The number of nodes}

  An object of class \code{\link[graph]{graphNEL-class}} is returned.
      
}

\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link{inducedGraph}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}

\keyword{graphs}
