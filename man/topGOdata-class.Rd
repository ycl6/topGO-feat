\name{topGOdata-class}

\docType{class}

\alias{topGOdata-class}

\alias{allGenes}
\alias{attrInTerm}
\alias{countGenesInTerm}
\alias{description<-}
\alias{description}
\alias{feasible<-}
\alias{feasible}
\alias{geneSelectionFun<-}
\alias{geneSelectionFun}
\alias{genes}
\alias{genesInTerm}

\alias{graph<-}
\alias{graph}
\alias{numGenes}
\alias{ontology<-}
\alias{ontology}
\alias{sigGenes}
\alias{numSigGenes}
\alias{termStat}
\alias{updateGenes}
\alias{updateTerm<-}
\alias{usedGO}
\alias{expressionMatrix}
\alias{phenotype}
\alias{expressionMatrix,topGOdata-method}
\alias{phenotype,topGOdata-method}

\alias{geneScore}
\alias{geneScore,topGOdata,missing-method}
\alias{geneScore,topGOdata,character-method}

\alias{scoresInTerm}
\alias{scoresInTerm,topGOdata,missing-method}
\alias{scoresInTerm,topGOdata,character-method}

\alias{show,topGOdata-method}
\alias{allGenes,topGOdata-method}
\alias{attrInTerm,topGOdata,character,character-method}
\alias{attrInTerm,topGOdata,character,missing-method}
\alias{countGenesInTerm,topGOdata,character-method}
\alias{countGenesInTerm,topGOdata,missing-method}
\alias{description<-,topGOdata,ANY-method}
\alias{description,topGOdata-method}
\alias{feasible<-,topGOdata-method}
\alias{feasible,topGOdata-method}
\alias{geneScore,topGOdata-method}
\alias{geneSelectionFun<-,topGOdata-method}
\alias{geneSelectionFun,topGOdata-method}
\alias{genes,topGOdata-method}
\alias{genesInTerm,topGOdata,character-method}
\alias{genesInTerm,topGOdata,missing-method}
\alias{graph<-,topGOdata-method}
\alias{graph,topGOdata-method}
\alias{initialize,topGOdata-method}
\alias{numGenes,topGOdata-method}
\alias{ontology<-,topGOdata-method}
\alias{ontology,topGOdata-method}
\alias{print,topGOdata-method}
\alias{sigGenes,topGOdata-method}
\alias{numSigGenes,topGOdata-method}
\alias{termStat,topGOdata,character-method}
\alias{termStat,topGOdata,missing-method}
\alias{updateGenes,topGOdata,numeric,function-method}
\alias{updateGenes,topGOdata,factor,missing-method}
\alias{updateTerm<-,topGOdata,character-method}
\alias{usedGO,topGOdata-method}

\title{Class "topGOdata"}
\description{
  TODO: The node attributes are environments containing the genes/probes
  annotated to the respective node

  If genes is a numeric vector than this should represent the gene's
  score. If it is factor it should discriminate the genes in
  interesting genes and the rest
  
  TODO: it will be a good idea to replace the allGenes and allScore with an
  ExpressionSet class. In this way we can use tests like global test, globalAncova....
  -- ALL variables starting with . are just for internal class usage (private)
  
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("topGOdata", ontology, allGenes, geneSelectionFun, description, annotationFun, ...)}.
  ~~ describe objects here ~~ 
}

\section{Slots}{
  \describe{
    \item{\code{description}:}{Object of class \code{"character"} ~~ }
    \item{\code{ontology}:}{Object of class \code{"character"} ~~ }
    \item{\code{allGenes}:}{Object of class \code{"character"} ~~ }
    \item{\code{allScores}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{geneSelectionFun}:}{Object of class \code{"function"} ~~ }
    \item{\code{feasible}:}{Object of class \code{"logical"} ~~ }
    \item{\code{nodeSize}:}{Object of class \code{"integer"} ~~ }
    \item{\code{graph}:}{Object of class \code{"graphNEL"} ~~ }
    \item{\code{expressionMatrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{phenotype}:}{Object of class \code{"factor"} ~~ }
    
  }
}

\section{Methods}{
  \describe{
    \item{allGenes}{\code{signature(object = "topGOdata")}: ... }
    \item{attrInTerm}{\code{signature(object = "topGOdata", attr = "character", whichGO = "character")}: ... }
    \item{attrInTerm}{\code{signature(object = "topGOdata", attr = "character", whichGO = "missing")}: ... }
    \item{countGenesInTerm}{\code{signature(object = "topGOdata", whichGO = "character")}: ... }
    \item{countGenesInTerm}{\code{signature(object = "topGOdata", whichGO = "missing")}: ... }
    \item{description<-}{\code{signature(object = "topGOdata")}: ... }
    \item{description}{\code{signature(object = "topGOdata")}: ... }
    \item{feasible<-}{\code{signature(object = "topGOdata")}: ... }
    \item{feasible}{\code{signature(object = "topGOdata")}: ... }
    \item{geneScore}{\code{signature(object = "topGOdata")}: ... }
    \item{geneSelectionFun<-}{\code{signature(object = "topGOdata")}: ... }
    \item{geneSelectionFun}{\code{signature(object = "topGOdata")}: ... }

    \item{genes}{\code{signature(object = "topGOdata")}: A method for
      obtaining the list of genes, as a characther vector, which will be
      used in the further analysis.}
    \item{numGenes}{\code{signature(object = "topGOdata")}: A method for
      obtaining the number of genes, which will be used in the further
      analysis. It has the same effect as: \code{lenght(genes(object))}.}

    \item{sigGenes}{\code{signature(object = "topGOdata")}:  A method for
      obtaining the list of significant genes, as a charachter vector.}
        
    \item{genesInTerm}{\code{signature(object = "topGOdata", whichGO = "character")}: ... }
    \item{genesInTerm}{\code{signature(object = "topGOdata", whichGO = "missing")}: ... }
    \item{getSigGroups}{\code{signature(object = "topGOdata", test.stat = "classicCount")}: ... }
    \item{getSigGroups}{\code{signature(object = "topGOdata", test.stat = "classicScore")}: ... }
    \item{graph<-}{\code{signature(object = "topGOdata")}: ... }
    \item{graph}{\code{signature(object = "topGOdata")}: ... }
    \item{initialize}{\code{signature(.Object = "topGOdata")}: ... }

    \item{ontology<-}{\code{signature(object = "topGOdata")}: ... }
    \item{ontology}{\code{signature(object = "topGOdata")}: ... }

    \item{termStat}{\code{signature(object = "topGOdata", whichGO = "character")}: ... }
    \item{termStat}{\code{signature(object = "topGOdata", whichGO = "missing")}: ... }
    \item{updateGenes}{\code{signature(object = "topGOdata", geneList = "numeric", geneSelFun = "function")}: ... }
    \item{updateGenes}{\code{signature(object = "topGOdata", geneList = "factor", geneSelFun = "missing")}: ... }
    \item{updateTerm<-}{\code{signature(object = "topGOdata", attr = "character")}: ... }
    \item{usedGO}{\code{signature(object = "topGOdata")}: ... }
  }
}

\author{Adrian Alexa}

\seealso{
  \code{\link{buildLevels}},
  \code{\link{annFUN}}
}

\examples{
## load the dataset 
data(geneList)
library(package = affyLib, character.only = TRUE)

## the distribution of the adjusted p-values
hist(geneList, 100)

## how many differentially expressed genes are:
sum(topDiffGenes(geneList))

## build the topGOdata class 
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSel = topDiffGenes,
              description = "GO analysis of ALL data: Differential Expression between 
                             B-cell and T-cell",
              annot = annFUN.db,
              affyLib = affyLib)

## display the GOdata object
GOdata

##########################################################
## Examples on how to use the methods
##########################################################

## description of the experiment
description(GOdata)

## obtain the genes that will be used in the analysis
a <- genes(GOdata)
str(a)
numGenes(GOdata)

## obtain the score (p-value) of the genes
selGenes <- names(geneList)[sample(1:length(geneList), 10)]
gs <- geneScore(GOdata, whichGenes = selGenes)
print(gs)

## if we want an unnamed vector containing all the feasible genes
gs <- geneScore(GOdata, use.names = FALSE)
str(gs)

## the list of significant genes
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)

## to update the gene list 
.geneList <- geneScore(GOdata, use.names = TRUE)
GOdata ## more available genes
GOdata <- updateGenes(GOdata, .geneList, topDiffGenes)
GOdata ## the available genes are now the feasible genes

## the available GO terms (all the nodes in the graph)
go <- usedGO(GOdata)
length(go)

## to list the genes annotated to a set of specified GO terms
sel.terms <- sample(go, 10)
ann.genes <- genesInTerm(GOdata, sel.terms)
str(ann.genes)

## the score for these genes
ann.score <- scoresInTerm(GOdata, sel.terms)
str(ann.score)

## to see the number of annotated genes
num.ann.genes <- countGenesInTerm(GOdata)
str(num.ann.genes)

## to summarise the statistics
termStat(GOdata, sel.terms)
}

\keyword{graphs}
\keyword{classes}
