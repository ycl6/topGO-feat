\name{getSigGroups}

\alias{getSigGroups}
\alias{runTest}

\alias{getSigGroups-methods}
\alias{getSigGroups,topGOdata,classicCount-method}
\alias{getSigGroups,topGOdata,classicScore-method}
\alias{getSigGroups,topGOdata,classicExpr-method}
\alias{getSigGroups,topGOdata,leaCount-method}
\alias{getSigGroups,topGOdata,leaScore-method}
\alias{getSigGroups,topGOdata,leaExpr-method}
\alias{getSigGroups,topGOdata,elimCount-method}
\alias{getSigGroups,topGOdata,elimScore-method}
\alias{getSigGroups,topGOdata,elimExpr-method}
\alias{getSigGroups,topGOdata,weightCount-method}
\alias{getSigGroups,topGOdata,weight01Count-method}
\alias{getSigGroups,topGOdata,weight01Score-method}
\alias{getSigGroups,topGOdata,weight01Expr-method}
\alias{getSigGroups,topGOdata,parentChild-method}
\alias{getSigGroups,topGOdata,pC-method}

\alias{runTest,topGOdata,character,character-method}
\alias{runTest,topGOdata,missing,character-method}


\alias{whichAlgorithms}
\alias{whichTests}


\title{Interfaces for running the enrichment tests}

\description{
  These function are used for dispatching the specific algorithm for
  a given \code{topGOdata} object and a test statistic.
}

\usage{
getSigGroups(object, test.stat, ...)
runTest(object, algorithm, statistic, ...)
whichAlgorithms()
whichTests()
}


\arguments{
  \item{object}{An object of class \code{topGOdata}
    This object contains all data necessary for runnig the test.}

  \item{test.stat}{An object of class \code{groupStats}. This object
    defines the test statistic.}
  
  \item{algorithm}{Character string specifing which algorithm to use.} 
  
  \item{statistic}{Character string specifing which test to use.} 
  
  \item{\dots}{Other parameters. In the case of \code{runTest} they are
    used for defining the test statistic}
}


\details{
  
  The \code{runTest} function can be used only with a predefined set of
  test statistics and algorithms. The algorithms and the statistical
  tests which are accessible via the \code{runTest} function are shown
  by the \code{whichAlgorithms()} and \code{whichTests()} functions.

  The runTest function is a warping of the \code{getSigGroups} and the
  initialisation of a \code{groupStats} object functions.
  
  ...
}


\value{
  An object of class \code{topGOresult}.
}


\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link{groupStats-class}},
  \code{\link{topGOresult-class}}
}



\examples{

## load a sample topGOdata object
data(GOdata)
GOdata

##############################
## getSigGroups interface
##############################

## define a test statistic
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
## perform the test
resultFis <- getSigGroups(GOdata, test.stat)
resultFis



##############################
## runTest interface
##############################

## Enrichment analysis by using the "classic" method and Fisher's exact test
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFis

## weight01 is the default algorithm 
weight01.fisher <- runTest(GOdata, statistic = "fisher")
weight01.fisher


## not all combinations are possible!
# weight.ks <- runTest(GOdata, algorithm = "weight", statistic = "t")
}

\keyword{methods}
