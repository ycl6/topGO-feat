\name{classicScore-class}
\docType{class}
\alias{classicScore-class}

\alias{scoreOrder}
\alias{membersScore}
\alias{rankMembers}
\alias{score<-}

\alias{allScore}
\alias{allScore,classicScore,missing-method}
\alias{allScore,classicScore,logical-method}

\alias{scoreOrder,classicScore-method}
\alias{initialize,classicScore-method}
\alias{membersScore,classicScore-method}
\alias{rankMembers,classicScore-method}
\alias{score<-,classicScore-method}


\alias{GOKSTest,classicScore-method}
\alias{GOtTest,classicScore-method}
\alias{GOSumTest,classicScore-method}
\alias{GOKSTiesTest,classicScore-method}



\title{Class "classicScore"}
\description{A class that extends the virtual class "groupStats" by adding 
  a slot representing the score of each gene. It is used for tests like
  Kolmogorov-Smirnov test.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("classicScore", testStatistic, name, allMembers, groupMembers, score, decreasing)}.
}

\section{Slots}{
	 \describe{
    \item{\code{score}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
    \item{\code{scoreOrder}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatPar}:}{Object of class \code{"ANY"} ~~ }
  }
}
\section{Extends}{
Class \code{"groupStats"}, directly.
}

\section{Methods}{
  \describe{
    \item{allScore}{Method to obtain the score of all members.}
    \item{scoreOrder}{Returns TRUE if the score should be ordered
      increasing, FALSE otherwise.}
    \item{membersScore}{\code{signature(object = "classicScore")}: ... }
    \item{rankMembers}{\code{signature(object = "classicScore")}: ... }
    \item{score<-}{\code{signature(object = "classicScore")}: ... }
	 }
}

\author{Adrian Alexa}

\seealso{
  \code{\link{classicCount-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}

\examples{
## define the type of test you want to use
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
}
\keyword{classes}
