\name{elimScore-class}
\docType{class}

\alias{leaScore-class}
\alias{elimScore-class}
\alias{weight01Score-class}

\alias{allMembers,elimScore-method}

\alias{allScore,elimScore,logical-method}
\alias{allScore,elimScore,missing-method}

\alias{cutOff,elimScore-method}
\alias{cutOff<-,elimScore-method}

\alias{elim,elimScore-method}
\alias{elim<-,elimScore-method}

\alias{depth,leaScore-method}
\alias{depth<-,leaScore-method}

\alias{alternative,elimScore-method}
\alias{membersScore,elimScore-method}
\alias{members,elimScore-method}
\alias{numMembers,elimScore-method}
\alias{numAllMembers,elimScore-method}
\alias{rankMembers,elimScore-method}
\alias{score<-,elimScore-method}


\alias{initialize,leaScore-method}
\alias{initialize,elimScore-method}
\alias{initialize,weight01Score-method}

\alias{allMembers,weight01Score-method}
\alias{allScore,weight01Score,missing-method}
\alias{allScore,weight01Score,logical-method}
\alias{elim<-,weight01Score-method}
\alias{elim,weight01Score-method}
\alias{members,weight01Score,missing-method}
\alias{membersScore,weight01Score-method}
\alias{numAllMembers,weight01Score-method}
\alias{numMembers,weight01Score-method}
\alias{rankMembers,weight01Score-method}




\title{Classes "elimScore" and "weight01Score" }
\description{Classes that extend the "classicScore" class by adding 
  a slot representing the members that need to be removed.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("elimScore", testStatistic, name, allMembers, groupMembers, score, alternative, elim, cutOff, ...)}.
	 ~~ describe objects here ~~ 
}


\section{Slots}{
	 \describe{
    \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{score}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.alternative}:}{Object of class \code{"logical"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
    \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{classicScore}"}, directly.
Class \code{"\linkS4class{groupStats}"}, by class "classicScore", distance 2.
}
\section{Methods}{
No methods defined with class "elimScore" in the signature.
}

\details{
  TODO:   
}

\author{Adrian Alexa}

\seealso{
  \code{\link{classicScore-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
}

\keyword{classes}
