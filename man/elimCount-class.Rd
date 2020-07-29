\name{elimCount-class}
\docType{class}

\alias{leaCount-class}
\alias{elimCount-class}
\alias{weight01Count-class}

\alias{cutOff}
\alias{cutOff,elimCount-method}
\alias{cutOff<-}
\alias{cutOff<-,elimCount-method}

\alias{elim}
\alias{elim,elimCount-method}
\alias{elim,weight01Count-method}
\alias{elim<-}
\alias{elim<-,elimCount-method}
\alias{elim<-,weight01Count-method}


\alias{depth}
\alias{depth,leaCount-method}
\alias{depth<-}
\alias{depth<-,leaCount-method}

\alias{initialize,leaCount-method}
\alias{initialize,elimCount-method}
\alias{initialize,weight01Count-method}

\alias{contTable,elimCount-method}
\alias{numSigAll,elimCount-method}
\alias{numSigMembers,elimCount-method}
\alias{numMembers,elimCount-method}
\alias{numAllMembers,elimCount-method}
\alias{sigAllMembers,elimCount-method}
\alias{sigMembers<-,elimCount-method}
\alias{sigMembers,elimCount-method}


\alias{numAllMembers,weight01Count-method}
\alias{numMembers,weight01Count-method}
\alias{numSigAll,weight01Count-method}
\alias{numSigMembers,weight01Count-method}
\alias{sigAllMembers,weight01Count-method}
\alias{sigMembers,weight01Count-method}

\alias{GOFisherTest,elimCount-method}


  
\title{Classes "elimCount" and "weight01Count"}
\description{Classes that extend the "classicCount" class by adding 
  a slot representing the members that need to be removed.}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("elimCount", testStatistic, name, allMembers, groupMembers, sigMembers, elim, cutOff, ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
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
  No methods defined with class "elimCount" in the signature.
}

\author{Adrian Alexa}

\seealso{
  \code{\link{classicScore-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}

\keyword{classes}
