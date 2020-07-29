\name{groupStats-class}
\docType{class}
\alias{groupStats-class}

\alias{allMembers<-}
\alias{allMembers}
\alias{members<-}
\alias{members}
\alias{Name<-}
\alias{Name}
\alias{numAllMembers}
\alias{numMembers}
\alias{testStatistic}
\alias{testStatPar}
\alias{updateGroup}



\alias{allMembers<-,groupStats-method}
\alias{allMembers,groupStats-method}
\alias{initialize,groupStats-method}
\alias{members<-,groupStats-method}
\alias{members,groupStats,missing-method}
\alias{Name<-,groupStats-method}
\alias{Name,groupStats-method}
\alias{numAllMembers,groupStats-method}
\alias{numMembers,groupStats-method}

\alias{runTest,groupStats-method}
\alias{runTest,groupStats,missing,missing-method}

\alias{testStatistic,groupStats-method}
\alias{testStatPar,groupStats-method}
\alias{updateGroup,groupStats,character,character-method}


\title{ Class "groupStats"}
\description{A virtual class containing basic gene set information:
  the gene universe, the member of the current group, the test statistic
  defined for this group, etc. 
}

\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
    \item{\code{testStatPar}:}{Object of class \code{"ANY"} ~~ }      
  }
}

\section{Methods}{
  \describe{
    \item{allMembers<-}{\code{signature(object = "groupStats")}: ... }
    \item{allMembers}{\code{signature(object = "groupStats")}: ... }
    \item{initialize}{\code{signature(.Object = "groupStats")}: ... }
    \item{members<-}{\code{signature(object = "groupStats")}: ... }
    \item{members}{\code{signature(object = "groupStats")}: ... }
    \item{Name<-}{\code{signature(object = "groupStats")}: ... }
    \item{Name}{\code{signature(object = "groupStats")}: ... }
    \item{numAllMembers}{\code{signature(object = "groupStats")}: ... }
    \item{numMembers}{\code{signature(object = "groupStats")}: ... }
    \item{runTest}{\code{signature(object = "groupStats")}: ... }
    \item{testStatistic}{\code{signature(object = "groupStats")}: ... }
	 }
}

\author{Adrian Alexa}

\seealso{
  \code{\link{classicCount-class}},
  \code{\link{getSigGroups-methods}}
}

\keyword{classes}
