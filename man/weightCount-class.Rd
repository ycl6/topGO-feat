\name{weightCount-class}
\docType{class}

\alias{weightCount-class}

\alias{getSigRatio}
\alias{sigRatio}
\alias{penalise}
\alias{roundFun}
\alias{significant}
\alias{sigRatio<-}

\alias{Weights}
\alias{Weights<-}
\alias{Weights,weightCount-method}
\alias{Weights,weightCount,missing-method}
\alias{Weights,weightCount,logical-method}
\alias{Weights<-,weightCount-method}

\alias{sigRatio,weightCount-method}
\alias{sigRatio<-,weightCount-method}

\alias{penalise,weightCount,numeric,numeric-method}

\alias{roundFun,weightCount-method}
\alias{significant,weightCount-method}
\alias{Name,weightCount-method}
\alias{allMembers,weightCount-method}
\alias{members,weightCount-method}
\alias{testStatistic,weightCount-method}
\alias{testStatPar,weightCount-method}
\alias{updateGroup,weightCount,character,character-method}

\alias{getSigRatio,weightCount-method}

\alias{initialize,weightCount-method}
\alias{numSigAll,weightCount-method}
\alias{numMembers,weightCount-method}
\alias{numSigMembers,weightCount-method}
\alias{numAllMembers,weightCount-method}


\title{Class "weightCount"}


\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}


\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("weightCount", testStatistic, name, allMembers, groupMembers, sigMembers, weights, sigRatio, penalise, ...)}.
}


\section{Slots}{
	 \describe{
    \item{\code{weights}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sigRatio}:}{Object of class \code{"function"} ~~ }
    \item{\code{penalise}:}{Object of class \code{"function"} ~~ }
    \item{\code{roundFun}:}{Object of class \code{"function"} ~~ }
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
  No methods defined with class "weightCount" in the signature.
}


\details{
  TODO: Some details here.....
  }

\author{Adrian Alexa}

\seealso{
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}


\keyword{classes}
