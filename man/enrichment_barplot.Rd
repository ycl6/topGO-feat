\name{enrichment_barplot}
\alias{enrichment_barplot}
\title{Create a barplot using ggplot2 to show enrichment results}
\usage{
enrichment_barplot(object, result, showTerms = 10, numChar = 40, 
		   orderBy = "Scores", y = "Count", 
		   xlab = NULL, ylab = NULL, title = NULL)
}
\arguments{
\item{object}{(Required). A \code{topGOdata} object.}

\item{result}{(Required). A \code{topGOresult} object.}

\item{showTerms}{(Optional). A single integer or a character vector of GO ID.
Default is \code{10}. Indicates the number of GO terms to show.}

\item{numChar}{(Optional). A single integer. Default is \code{40}.
Indicates the number characters to keep in the GO term description}

\item{orderBy}{(Optional). A character string. Default is \code{"Scores"}.
Indicates how to order the test results before subsetting to taking top \code{N} terms.
It can be either \code{"Scores"} (i.e. \code{p}-values) or \code{"Ratio"}}

\item{y}{(Optional). A character string. Default is \code{"Count"}.
Indicates the variable that should be mapped to the y-axis.
It can be either \code{"Count"} or \code{"Ratio"}}

\item{xlab}{(Optional). A character string. Default is \code{NULL}.
Indicates the x-axis label.}

\item{ylab}{(Optional). A character string. Default is \code{NULL}.
Indicates the y-axis label.}

\item{title}{(Optional). A character string. Default is \code{NULL}
Indicates the main title for the graphic.}
}
\value{
A \code{\link{ggplot}}2 plot object
}
\description{
This function takes a \code{topGOdata} object and a \code{topGOresult} object
to display the enrichment statistics of top \code{N} GO terms.
It wraps \code{ggplot2} plotting, and returns a \code{ggplot2} graphic object.
}
\examples{
data(GOdata)
data(results.tGO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
enrichment_barplot(GOdata, resultFisher, showTerms = 10, numChar = 40,
                   orderBy = "Scores", y = "Count")
}
\seealso{
\code{\link[ggplot2]{ggplot}}
}

