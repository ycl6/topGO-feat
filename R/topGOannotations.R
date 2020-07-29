########################################################
##
## functions related to annotations ...
## 
########################################################

readMappings <- function(file, sep = "\t", IDsep = ",") {
  a <- read.delim(file = file, header = FALSE,
                  quote = "", sep = sep, colClasses = "character")
  
  ## a bit of preprocesing to get it in a nicer form
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1]) ## trim the spaces
  
  ## split the IDs
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]])))
}


## function annFUN.org() to work with the "org.XX.eg" annotations
annFUN.org <- function(whichOnto, feasibleGenes = NULL, mapping, ID = "entrez") {

  tableName <- c("genes", "accessions", "alias", "ensembl",
                 "gene_info", "gene_info", "unigene")
  keyName <- c("gene_id", "accessions", "alias_symbol", "ensembl_id",
               "symbol", "gene_name", "unigene_id")
  names(tableName) <- names(keyName) <- c("entrez", "genbank", "alias", "ensembl",
                                          "symbol", "genename", "unigene")
  
  
  ## we add the .db ending if needed 
  mapping <- paste(sub(".db$", "", mapping), ".db", sep = "")
  require(mapping, character.only = TRUE) || stop(paste("package", mapping, "is required", sep = " "))
  mapping <- sub(".db$", "", mapping)
  
  geneID <- keyName[tolower(ID)]
  .sql <- paste("SELECT DISTINCT ", geneID, ", go_id FROM ", tableName[tolower(ID)],
                " INNER JOIN ", paste("go", tolower(whichOnto), sep = "_"),
                " USING(_id)", sep = "")
  retVal <- dbGetQuery(get(paste(mapping, "dbconn", sep = "_"))(), .sql)
  
  ## restric to the set of feasibleGenes
  if(!is.null(feasibleGenes))
    retVal <- retVal[retVal[[geneID]] %in% feasibleGenes, ]
  
  ## split the table into a named list of GOs
  return(split(retVal[[geneID]], retVal[["go_id"]]))
}


annFUN.db <- function(whichOnto, feasibleGenes = NULL, affyLib) {

  ## we add the .db ending if needed 
  affyLib <- paste(sub(".db$", "", affyLib), ".db", sep = "")
  require(affyLib, character.only = TRUE) || stop(paste("package", affyLib, "is required", sep = " "))
  affyLib <- sub(".db$", "", affyLib)

  orgFile <- get(paste(get(paste(affyLib, "ORGPKG", sep = "")), "_dbfile", sep = ""))
  
  try(dbGetQuery(get(paste(affyLib, "dbconn", sep = "_"))(),
                 paste("ATTACH '", orgFile(), "' as org;", sep ="")),
      silent = TRUE)
  
  .sql <- paste("SELECT DISTINCT probe_id, go_id FROM probes INNER JOIN ",
                "(SELECT * FROM org.genes INNER JOIN org.go_",
                tolower(whichOnto)," USING('_id')) USING('gene_id');", sep = "")

  retVal <- dbGetQuery(get(paste(affyLib, "dbconn", sep = "_"))(), .sql)

  ## restric to the set of feasibleGenes
  if(!is.null(feasibleGenes))
    retVal <- retVal[retVal[["probe_id"]] %in% feasibleGenes, ]
  
  ## split the table into a named list of GOs
  return(split(retVal[["probe_id"]], retVal[["go_id"]]))
}



annFUN <- function(whichOnto, feasibleGenes = NULL, affyLib) {
    
  require(affyLib, character.only = TRUE) || stop(paste('package', affyLib, 'is required', sep = " "))
  mapping <- get(paste(affyLib, 'GO2PROBE', sep = ''))

  if(is.null(feasibleGenes))
    feasibleGenes <- ls(get(paste(affyLib, 'ACCNUM', sep = '')))
      
  ontoGO <- get(paste('GO', whichOnto, "Term", sep = ''))
  goodGO <- intersect(ls(ontoGO), ls(mapping))

  GOtoAffy <- lapply(mget(goodGO, envir = mapping, ifnotfound = NA),
                     intersect, feasibleGenes)
  
  emptyTerms <- sapply(GOtoAffy, length) == 0

  return(GOtoAffy[!emptyTerms])
}


## the annotation function
annFUN.gene2GO <- function(whichOnto, feasibleGenes = NULL, gene2GO) {

  ## GO terms annotated to the specified ontology 
  ontoGO <- get(paste("GO", whichOnto, "Term", sep = ""))

  ## Restrict the mappings to the feasibleGenes set  
  if(!is.null(feasibleGenes))
    gene2GO <- gene2GO[intersect(names(gene2GO), feasibleGenes)]
  
  ## Throw-up the genes which have no annotation
  if(any(is.na(gene2GO)))
    gene2GO <- gene2GO[!is.na(gene2GO)]
  
  gene2GO <- gene2GO[sapply(gene2GO, length) > 0]
  
  ## Get all the GO and keep a one-to-one mapping with the genes
  allGO <- unlist(gene2GO, use.names = FALSE)
  geneID <- rep(names(gene2GO), sapply(gene2GO, length))

  goodGO <- allGO %in% ls(ontoGO)
  return(split(geneID[goodGO], allGO[goodGO]))
}


## the annotation function
annFUN.GO2genes <- function(whichOnto, feasibleGenes = NULL, GO2genes) {
  
  ## GO terms annotated to the specified ontology 
  ontoGO <- get(paste("GO", whichOnto, "Term", sep = ""))

  ## Select only the GO's form the specified ontology
  GO2genes <- GO2genes[intersect(ls(ontoGO), names(GO2genes))]

  ## Restrict the mappings to the feasibleGenes set  
  if(!is.null(feasibleGenes))
    GO2genes <- lapply(GO2genes, intersect, feasibleGenes)

  return(GO2genes[sapply(GO2genes, length) > 0])
}


## annotation function to read the mappings from a file
annFUN.file <- function(whichOnto, feasibleGenes = NULL, file, ...) {

  ## read the mappings from the file
  mapping <- readMappings(file = file, ...)

  ## we check which direction the mapping is ...
  GOs <- ls(get(paste("GO", whichOnto, "Term", sep = "")))
  if(any(GOs %in% names(mapping)))
    return(annFUN.GO2genes(whichOnto, feasibleGenes, mapping))
  
  return(annFUN.gene2GO(whichOnto, feasibleGenes, mapping))
}  


## function returning all genes that can be used for analysis 
feasibleGenes.db <- function(affyLib, whichOnto) {

  affyLib <- sub(".db$", "", affyLib)
  mapping <- get(paste(affyLib, 'GO2PROBE', sep = ''))
  
  ontoGO <- get(paste('GO', whichOnto, "Term", sep = ''))
  goodGO <- intersect(ls(ontoGO), ls(mapping))

  return(unique(unlist(mget(goodGO, envir = mapping, ifnotfound = NA))))
}

