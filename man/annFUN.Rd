\name{annFUN}

\alias{annFUN.db}
\alias{annFUN}
\alias{annFUN.GO2genes}
\alias{annFUN.gene2GO}
\alias{annFUN.file}
\alias{annFUN.org}

\alias{inverseList}
\alias{readMappings}


\title{Functions which map gene identifiers to GO terms}
\description{
  These functions are used to compile a list of GO terms such that each
  element in the list is a character vector containing all the gene
  identifiers that are mapped to the respective GO term.  
}


\usage{
annFUN.db(whichOnto, feasibleGenes = NULL, affyLib)
annFUN.org(whichOnto, feasibleGenes = NULL, mapping, ID = "entrez") 
annFUN(whichOnto, feasibleGenes = NULL, affyLib)
annFUN.gene2GO(whichOnto, feasibleGenes = NULL, gene2GO)
annFUN.GO2genes(whichOnto, feasibleGenes = NULL, GO2genes)
annFUN.file(whichOnto, feasibleGenes = NULL, file, ...)

readMappings(file, sep = "\t", IDsep = ",")
inverseList(l)
}


\arguments{
  \item{whichOnto}{character string specifying one of the three GO
    ontologies, namely: \code{"BP"}, \code{"MF"}, \code{"CC"}}

  \item{feasibleGenes}{character vector containing a subset of gene
    identifiers. Only these genes will be used to annotate GO
    terms. Default value is \code{NULL} which means that there are no
    genes filtered.}

  \item{affyLib}{character string containing the name of the
    Bioconductor annotaion package for a specific microarray chip.}

  \item{gene2GO}{named list of character vectors. The list names are
    genes identifiers. For each gene the character vector contains the
    GO identifiers it maps to. Only the most specific annotations are required.}
  
  \item{GO2genes}{named list of character vectors. The list names are
    GO identifiers. For each GO the character vector contains the
    genes identifiers which are mapped to it. Only the most specific
    annotations are required.}
  
  \item{mapping}{character string specifieng the name of the
    Bioconductor package containing the gene mappings for a
    specific organism. For example: \code{mapping = "org.Hs.eg.db"}.}
  
  \item{ID}{character string specifing the gene identifier to
    use. Currently only the following identifiers can be used:
    \code{c("entrez", "genbank", "alias", "ensembl", "symbol",
      "genename", "unigene")}}
  
  \item{file}{character string specifing the file containing the annotations.}

  \item{...}{other parameters}

  \item{sep}{the character used to separate the columns in the CSV file}
  
  \item{IDsep}{the character used to separate the annotated entities}

  \item{l}{a list containing mappings}
}

\details{
  All these function restrict the GO terms to the ones belonging
  to the specified ontology and to the genes listed in the
  \code{feasibleGenes} attribute (if not empty).
  
  The function \code{annFUN.db} uses the mappings provided
  in the Bioconductor annotation data packages. For example, if the
  Affymetrix \code{hgu133a} chip it is used, then the user should set
  \code{affyLib = "hgu133a.db"}.

  The functions \code{annFUN.gene2GO} and \code{annFUN.GO2genes} are
  used when the user provide his own annotations either as a gene-to-GOs
  mapping, either as a GO-to-genes mapping.
  
  The \code{annFUN.org} function is using the mappings from the
  "org.XX.XX" annotation packages. The function supports different gene
  identifiers.

  The \code{annFUN.file} function will read the annotationsof the type
  gene2GO or GO2genes from a text file.          
}


\value{
  A named(GO identifiers) list of character vectors. 
}


\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}}
}

\examples{

library(hgu133a.db)
set.seed(111)

## generate a gene list and the GO annotations
selGenes <- sample(ls(hgu133aGO), 50)
gene2GO <- lapply(mget(selGenes, envir = hgu133aGO), names)
gene2GO[sapply(gene2GO, is.null)] <- NA

## the annotation for the first three genes
gene2GO[1:3]

## inverting the annotations
G2g <- inverseList(gene2GO)

## inverting the annotations and selecting an ontology
go2genes <- annFUN.gene2GO(whichOnto = "CC", gene2GO = gene2GO)


## generate a GO list with the genes annotations
selGO <- sample(ls(hgu133aGO2PROBE), 30)
GO2gene <- lapply(mget(selGO, envir = hgu133aGO2PROBE), as.character)

GO2gene[1:3]

## select only the GO terms for a specific ontology
go2gene <- annFUN.GO2genes(whichOnto = "CC", GO2gene = GO2gene)


##################################################
## Using the org.XX.xx.db annotations
##################################################

## GO to Symbol mappings (only the BP ontology is used)
xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
head(xx)

\dontrun{

allGenes <- unique(unlist(xx))
myInterestedGenes <- sample(allGenes, 500)
geneList <- factor(as.integer(allGenes %in% myInterestedGenes))
names(geneList) <- allGenes

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              nodeSize = 5,
              annot = annFUN.org, 
              mapping = "org.Hs.eg.db",
              ID = "symbol") 
}

}

\keyword{misc}
