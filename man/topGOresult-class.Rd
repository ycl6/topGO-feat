\name{topGOresult-class}
\docType{class}

\alias{topGOresult-class}


\alias{initialize,topGOresult-method}

\alias{description,topGOresult-method}
\alias{description<-,topGOresult,ANY-method}

\alias{score}
\alias{score,topGOresult-method}
\alias{score<-,topGOresult-method}

\alias{testName}
\alias{testName<-}
\alias{testName,topGOresult-method}
\alias{testName<-,topGOresult-method}

\alias{algorithm}
\alias{algorithm<-}
\alias{algorithm<-,topGOresult-method}
\alias{algorithm,topGOresult-method}

\alias{show,topGOresult-method}
\alias{print,topGOresult-method}

\alias{geneData<-}
\alias{geneData}
\alias{geneData<-,topGOresult-method}
\alias{geneData,topGOresult-method}

\alias{combineResults}





\title{Class "topGOresult"}

\description{Class instance created by
  \code{\link{getSigGroups-methods}} or by \code{runTest}}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("topGOresult",
    description, score, testName, algorithm, geneData)}.
}

\section{Slots}{
  \describe{
    \item{\code{description}:}{character string containing a
      short description on how the object was build.}
    \item{\code{score}:}{named numerical vector containing the
      p-values or the scores of the tested GO terms.}
    \item{\code{testName}:}{character string containing the name
      of the test statistic used.}
    \item{\code{algorithm}:}{character string containing the name
      of the algorithm used.}
    \item{\code{geneData}:}{list containing summary statistics
      on the genes/gene universe/annotations.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{score}:}{method to access the \code{score} slot.}
    
    \item{\code{testName}:}{method to access the \code{testName} slot.}

    \item{\code{algorithm}:}{method to access the \code{algorithm} slot.}

    \item{\code{geneData}:}{method to access the \code{geneData} slot.}

    \item{\code{show}:}{method to print the object.}

    \item{\code{combineResults}:}{method to aggregate two or more
      topGOresult objects. \code{method = c("gmean", "mean", "median",
      "min", "max")} provides the way the object scores (which most of
      the time are p-values) are combined.}.
}}


\author{Adrian Alexa}

\seealso{
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}

\examples{

data(results.tGO)

s <- score(resultFisher)

go <- sort(names(s))
go.sub<- sample(go, 100)
go.mixed <- c(sample(go, 50), sample(ls(GOCCTerm), 20))
go.others <- sample(ls(GOCCTerm), 100)


str(go)
str(go.sub)
str(go.mixed)
str(go.others)

str(score(resultFisher, whichGO = go))
str(score(resultFisher, whichGO = go.sub))
str(score(resultFisher, whichGO = go.mixed))
str(score(resultFisher, whichGO = go.others))

avgResult <- combineResults(resultFisher, resultKS)
avgResult
combineResults(resultFisher, resultKS, method = "min")

}

\keyword{classes}
