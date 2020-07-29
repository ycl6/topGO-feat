\name{classicCount-class}
\docType{class}
\alias{classicCount-class}

\alias{contTable}
\alias{numSigAll}
\alias{numSigMembers}
\alias{sigAllMembers}
\alias{sigMembers<-}
\alias{sigMembers}

\alias{contTable,classicCount-method}
\alias{initialize,classicCount-method}
\alias{numSigAll,classicCount-method}
\alias{numSigMembers,classicCount-method}
\alias{sigAllMembers,classicCount-method}
\alias{sigMembers<-,classicCount-method}
\alias{sigMembers,classicCount-method}
\alias{GOFisherTest,classicCount-method}

\title{Class "classicCount"}
\description{This class that extends the virtual class "groupStats" by adding 
  a slot representing the significant members.}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("classicCount",
    testStatistic = "function",
    name = "character",
    allMembers = "character",
    groupMembers = "character",
    sigMembers = "character")}.
}

\section{Slots}{
	 \describe{
    \item{\code{significant}:}{Object of class \code{"integer"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
  }
}

\section{Extends}{
Class \code{"groupStats"}, directly.
}

\section{Methods}{
  \describe{
    \item{contTable}{\code{signature(object = "classicCount")}: ... }
    \item{initialize}{\code{signature(.Object = "classicCount")}: ... }
    \item{numSigAll}{\code{signature(object = "classicCount")}: ... }
    \item{numSigMembers}{\code{signature(object = "classicCount")}: ... }
    \item{sigAllMembers}{\code{signature(object = "classicCount")}: ... }
    \item{sigMembers<-}{\code{signature(object = "classicCount")}: ... }
    \item{sigMembers}{\code{signature(object = "classicCount")}: ... }
	 }
}

\details{
  This class is used for test statistic based on counts, like Fisher's
  exact test
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
