\name{classicExpr-class}
\docType{class}
\alias{classicExpr-class}

\alias{pType<-}
\alias{pType}
\alias{membersExpr}

\alias{allMembers<-,classicExpr-method}
\alias{emptyExpr,classicExpr-method}
\alias{GOglobalTest,classicExpr-method}
\alias{initialize,classicExpr-method}
\alias{membersExpr,classicExpr-method}
\alias{pType<-,classicExpr-method}
\alias{pType,classicExpr-method}

\title{Class "classicExpr"}
\description{This class that extends the virtual class "groupStats" by adding 
  two slots for accomodating gene expression data.}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("classicExpr", testStatistic, name, groupMembers, exprDat, pType, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{eData}:}{Object of class \code{"environment"} ~~ }
    \item{\code{pType}:}{Object of class \code{"factor"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
    \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{groupStats}"}, directly.
}
\section{Methods}{
  \describe{
    \item{allMembers<-}{\code{signature(object = "classicExpr")}: ... }
    \item{emptyExpr}{\code{signature(object = "classicExpr")}: ... }
    \item{getSigGroups}{\code{signature(object = "topGOdata", test.stat = "classicExpr")}: ... }
    \item{GOglobalTest}{\code{signature(object = "classicExpr")}: ... }
    \item{initialize}{\code{signature(.Object = "classicExpr")}: ... }
    \item{membersExpr}{\code{signature(object = "classicExpr")}: ... }
    \item{pType<-}{\code{signature(object = "classicExpr")}: ... }
    \item{pType}{\code{signature(object = "classicExpr")}: ... }
  }
}

\author{Adrian Alexa}

\seealso{
  \code{\link{classicScore-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}

\examples{
showClass("classicExpr")
}
\keyword{classes}
