\name{parentChild-class}
\docType{class}

\alias{parentChild-class}
\alias{pC-class}

\alias{allParents}
\alias{joinFun}

\alias{allMembers<-,parentChild-method}
\alias{allMembers,parentChild-method}
\alias{allParents,parentChild-method}
\alias{initialize,parentChild-method}
\alias{joinFun,parentChild-method}
\alias{numAllMembers,parentChild-method}
\alias{numSigAll,parentChild-method}
\alias{sigAllMembers,parentChild-method}
\alias{sigMembers<-,parentChild-method}
\alias{updateGroup,parentChild,missing,character-method}

\alias{allMembers<-,pC-method}
\alias{initialize,pC-method}
\alias{sigMembers<-,pC-method}
\alias{updateGroup,pC,missing,character-method}
\alias{updateGroup,pC,missing,missing-method}



\title{Classes "parentChild" and "pC"}
\description{Classes that extend the "classicCount" class by adding 
  support for the parent-child test.}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("parentChild", testStatistic, name, groupMembers, parents, sigMembers, joinFun, ...)}.
}

\section{Slots}{
	 \describe{
    \item{\code{splitIndex}:}{Object of class \code{"integer"} ~~ }
    \item{\code{joinFun}:}{Object of class \code{"character"} ~~ }
    \item{\code{significant}:}{Object of class \code{"integer"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
    \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{classicCount}"}, directly.
Class \code{"\linkS4class{groupStats}"}, by class "classicCount", distance 2.
}
\section{Methods}{
  \describe{
    \item{allMembers<-}{\code{signature(object = "parentChild")}: ... }
    \item{allMembers}{\code{signature(object = "parentChild")}: ... }
    \item{allParents}{\code{signature(object = "parentChild")}: ... }
    \item{getSigGroups}{\code{signature(object = "topGOdata", test.stat = "parentChild")}: ... }
    \item{initialize}{\code{signature(.Object = "parentChild")}: ... }
    \item{joinFun}{\code{signature(object = "parentChild")}: ... }
    \item{numAllMembers}{\code{signature(object = "parentChild")}: ... }
    \item{numSigAll}{\code{signature(object = "parentChild")}: ... }
    \item{sigAllMembers}{\code{signature(object = "parentChild")}: ... }
    \item{sigMembers<-}{\code{signature(object = "parentChild")}: ... }
    \item{updateGroup}{\code{signature(object = "parentChild", name = "missing", members = "character")}: ... }
	 }
}


\author{Adrian Alexa}

\seealso{
  \code{\link{classicCount-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}


\examples{
showClass("parentChild")
showClass("pC")
}

\keyword{classes}
