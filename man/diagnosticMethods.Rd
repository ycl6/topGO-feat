\name{dignostic-methods}

\alias{printGenes-methods}
\alias{printGenes}
\alias{printGenes,topGOdata,character,character-method}
\alias{printGenes,topGOdata,character,missing-method}

\alias{GenTable}
\alias{GenTable,topGOdata-method}

\alias{showGroupDensity}

\title{Diagnostic functions for topGOdata and topGOresult objects.}

\description{The \code{GenTable} function generates a summary of the
  results of the enrichment analysis.
  
  The \code{showGroupDensity} function plots the distributions of the
  gene' scores/ranks inside a GO term.

  The \code{printGenes} function shows a short summary of the top genes annotated
  to the specified GO terms.
}


\usage{
GenTable(object, ...)

showGroupDensity(object, whichGO, ranks = FALSE, rm.one = TRUE) 

printGenes(object, whichTerms, file, ...)
}

\arguments{
  \item{object}{an object of class \code{topGOdata}.}
  \item{whichGO}{the GO terms for which the plot should be generated.}
  \item{ranks}{if ranks should be used instead of scores.}
  \item{rm.one}{the p-values which are 1 are removed. }  
  \item{whichTerms}{
    character vector listing the GO terms for which
    the summary should be printed.
  }
  \item{file}{
    character string specifying the file in which the results should
    be printed.
  }
  \item{\dots}{
    %%GenTable(object, ..., orderBy = 1, ranksOf = 2, topNodes = 10, numChar = 40)
    Extra arguments for \code{GenTable} can be:
    \describe{
      \item{\dots}{
	one or more objects of class \code{topGOresult}.
      }
      \item{\code{orderBy}}{
	if more than one \code{topGOresult} object is given then
	\code{orderBy} gives the index of which scores will be
	used to order the resulting table. Can be an integer index
	or a character vector given the name of the \code{topGOresult}
	object.
      }
      \item{\code{ranksOf}}{
	same as \code{orderBy} argument except that this parameter shows
	the relative ranks of the specified result.
      }
      \item{\code{topNodes}}{
	the number of top GO terms to be included in the table.
      }
      \item{\code{numChar}}{
	the GO term definition will be truncated such that only
	the first \code{numChar} characters are shown.
      }
    }
    
    %% printGenes(object, whichTerms, chip, numChar = 100, simplify = TRUE, geneCufOff = 50, pvalCutOff)
    %% printGenes(object, whichTerms, file, oneFile = FALSE, ...)
    Extra arguments for \code{printGenes} can be:
    \describe{
      \item{\code{chip}}{
	character string containing the name of the Bioconductor
	annotation package for a microarray chip.
      }
      \item{\code{numChar}}{
	the gene description is trimmed such that it has
	\code{numChar} characters.
      }
      \item{\code{simplify}}{
	logical variable affecting how the results are returned.
      }
      \item{\code{geneCutOff}}{
	the maximal number of genes shown for each term.
      }
      \item{\code{pvalCutOff}}{
	only the genes with a p-value less than \code{pvalCutOff} are shown.
      }
      \item{\code{oneFile}}{
	if \code{TRUE} then a file for each GO term is generated.
      }
    }
  }
}
  
\details{

  \code{GenTable} is an easy to use function for summarising the most
  significant GO terms and the corresponding p-values. The function
  dispatches for \code{topGOdata} and \code{topGOresult} objects, and
  it can take an arbitrary number of the later, making comparison
  between various results easier.
  
  Note: One needs to type the complete attribute names (the exact name)
  of this function, like: \code{topNodes = 5}, \code{rankOf = "resultFis"}, etc. 
  This being the price paid for flexibility of specifying different
  number of \code{topGOdata} objects.
  

  
  The \code{showGroupDensity} function analyse the distribution of the
  gene-wise scores for a specified GO term.
  The function will show the distribution of the genes in a GO term
  compared with the complementary set, using a lattice plot.
  

  \code{printGenes} 
  The function will generate a table with all the probes annotated to
  the specified GO term. Various type of identifiers, the gene name and
  the gene-wise statistics are provided in the table. 
  
  One or more GO identifiers can be given to the function using the
  \code{whichTerms} argument. When more than one GO is specified, the
  function returns a list of \code{data.frames}, otherwise only one
  \code{data.frame} is returned.
  
  The function has a argument \code{file} which, when specified, will
  save the results into a file using the CSV format.

  For the moment the function will work only when the chip used has an
  annotation package available in Bioconductor. It will not work with
  other type of custom annotations.
}

\value{
  A data.frame or a list of data.fames.
}

\author{Adrian Alexa}

\seealso{
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}


\examples{

data(GOdata)


########################################
## GenTable
########################################

## load two topGOresult sample objects: resultFisher and resultKS
data(results.tGO)

## generate the result of Fisher's exact test
sig.tab <- GenTable(GOdata, Fis = resultFisher, topNodes = 20)

## results of both test
sig.tab <- GenTable(GOdata, resultFisher, resultKS, topNodes = 20)

## results of both test with specified names
sig.tab <- GenTable(GOdata, Fis = resultFisher, KS = resultKS, topNodes = 20)

## results of both test with specified names and specified ordering
sig.tab <- GenTable(GOdata, Fis = resultFisher, KS = resultKS, orderBy = "KS", 
		    ranksOf = "Fis", topNodes = 20)


########################################
## showGroupDensity
########################################

goID <- "GO:0006091"
print(showGroupDensity(GOdata, goID, ranks = TRUE))
print(showGroupDensity(GOdata, goID, ranks = FALSE, rm.one = FALSE))



########################################
## printGenes
########################################

\dontrun{
library(hgu95av2.db)
goID <- "GO:0006629"

gt <- printGenes(GOdata, whichTerms = goID, chip = "hgu95av2.db", numChar = 40)

goIDs <- c("GO:0006629", "GO:0007076")
gt <- printGenes(GOdata, whichTerms = goIDs, chip = "hgu95av2.db", pvalCutOff = 0.01)

gt[goIDs[1]]
}
}

\keyword{methods}

