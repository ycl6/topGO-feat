readMappings <- function(file, sep = "\t", IDsep = ",") {
  a <- read.delim(file = file, header = FALSE,
                  quote = "", sep = sep, colClasses = "character")
  
  ## a bit of preprocesing to get it in a nicer form
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1]) ## trim the spaces
  
  ## split the IDs
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]])))
}



gene2GO <- readMappings(file = "hs-ensembl2go.map")

set.seed(length(gene2GO))
index <- sample(length(gene2GO), 100)


str(head(gene2GO[index]))


## at this point we have the gene to GOs mapping and we can start building a topGOdata object.
## We need to generate a random list of interesting genes


## get the list of all gene ID's
geneNames <- names(gene2GO[index])

## select (or define) the list of interesting genes
myInterestedGenes <- sample(geneNames, 20)

str(geneNames)
str(myInterestedGenes)

## make a indicator vector showing which genes are interesting
geneList <- factor(as.integer(geneNames %in% myInterestedGenes))
names(geneList) <- geneNames

str(geneList)


library(topGO)

## build the topGOdata class
## there are three annotation functions available:
##      1. annFUN.db  -- used for bioconductor affy chips
##      2. annFUN.gene2GO  -- used when you have mappings from each gene to GOs
##      3. annFUN.GO2genes -- used when you have mappings from each GO to genes
##
help("annFUN")

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,  ## the new annotation function
              gene2GO = gene2GO[index])      ## the GO to gene ID's dataset






## write the gene2GO[index] into a file to use as an example in the vignette
a <- lapply(gene2GO[index], function(x) paste(x, collapse = ", "))
mat <- cbind(names(a), unlist(a, use.names = FALSE))
write.table(mat, file = "geneid2go.map", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

fn <- system.file("examples/geneid2go.map", package = "topGO")

