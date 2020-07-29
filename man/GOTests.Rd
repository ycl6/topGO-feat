\name{Gene set tests statistics}

\alias{GOFisherTest}
\alias{GOKSTest}
\alias{GOtTest}
\alias{GOglobalTest}
\alias{GOSumTest}
\alias{GOKSTiesTest}
\alias{permSumStats}
\alias{permSumStats.all}

\title{Gene set tests statistics}

\description{Methods which implement and run a group test statistic for
  a class inheriting from \code{groupStats} class. See Details section
  for a description of each method.
}

\usage{
GOFisherTest(object)
GOKSTest(object)
GOtTest(object)
GOglobalTest(object)
GOSumTest(object)
GOKSTiesTest(object)
}

\arguments{
  \item{object}{An object of class \code{groupStats} or decedent class.}
}

\details{
  GOFisherTest: implements Fischer's exact test (based on contingency
  table) for \code{groupStats} objects dealing with "counts".
  
  GOKSTest: implements the Kolmogorov-Smirnov test for \code{groupStats}
  objects dealing with gene "scores". This test uses the \code{ks.test}
  function and does not implement the running-sum-statistic test based
  on permutations. 

  GOtTest: implements the t-test for \code{groupStats} objects dealing
  with gene "scores". It should be used when the gene scores are
  t-statistics or any other score following a normal distribution.
    
  GOglobalTest: implement Goeman's globaltest.
}

\value{
  All these methods return the p-value computed by the respective test statistic.
}

\author{Adrian Alexa}

\seealso{
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}


\keyword{misc}
