\name{elimExpr-class}
\docType{class}

\alias{leaExpr-class}
\alias{elimExpr-class}
\alias{weight01Expr-class}

\alias{depth,leaExpr-method}
\alias{depth<-,leaExpr-method}

\alias{cutOff<-,elimExpr-method}
\alias{cutOff,elimExpr-method}

\alias{initialize,leaExpr-method}
\alias{initialize,elimExpr-method}
\alias{initialize,weight01Expr-method}


\alias{allMembers,weight01Expr-method}
\alias{elim<-,weight01Expr-method}
\alias{elim,weight01Expr-method}
\alias{members,weight01Expr,missing-method}
\alias{numAllMembers,weight01Expr-method}
\alias{numMembers,weight01Expr-method}


\title{Class "elimExpr"}
\description{Classes that extend the "classicExpr" class by adding 
  a slot representing the members that need to be removed.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("elimExpr", testStatistic, name, groupMembers, exprDat, pType, elim, cutOff, ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
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
Class \code{"\linkS4class{weight01Expr}"}, directly.
Class \code{"\linkS4class{classicExpr}"}, by class "weight01Expr", distance 2.
Class \code{"\linkS4class{groupStats}"}, by class "weight01Expr", distance 3.
}
\section{Methods}{
  \describe{
    \item{cutOff<-}{\code{signature(object = "elimExpr")}: ... }
    \item{cutOff}{\code{signature(object = "elimExpr")}: ... }
    \item{getSigGroups}{\code{signature(object = "topGOdata", test.stat = "elimExpr")}: ... }
    \item{initialize}{\code{signature(.Object = "elimExpr")}: ... }
	 }
}


\details{
  TODO: Some datails here.....
}

\author{Adrian Alexa}

\seealso{
  \code{\link{classicScore-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}

\examples{
showClass("elimExpr")
}
\keyword{classes}
