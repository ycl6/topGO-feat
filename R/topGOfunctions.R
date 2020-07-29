#################### misc functions ####################
##
##
##
##
##
########################################################



## Function that split GOTERM in different ontologies.
## Every new environment contain only the terms from one
## of the ontologies 'BP', 'CC', 'MF'
groupGOTerms <- function(where) {

  if(missing(where))
    where <- .GlobalEnv
  where <- as.environment(where)

  sql <- "SELECT go_id FROM go_term WHERE ontology IN"
  for(onto in c("BP", "MF", "CC")) {
    xx <- dbGetQuery(GO_dbconn(), paste(sql, "('", onto, "');", sep = ""))$go_id
    e <- new.env(hash = T, parent = emptyenv())
    multiassign(xx, value = rep(TRUE, length(xx)), envir = e)
    assign(paste("GO", onto, "Term", sep = ""), e, envir = where)
  }

  message("\ngroupGOTerms: \tGOBPTerm, GOMFTerm, GOCCTerm environments built.")
}


## given a list of vectors this function is returning the reverse list
## given a mapping form genes to GO terms as a list, compute which are
## the genes mapped to each GO.
inverseList <- function(l) {
  rId <- unlist(l, use.names = FALSE)
  lId <- rep(names(l), sapply(l, length))

  return(split(lId, rId))
}





######################################################################
## function to print all the genes annotated to the specified GOs

if(!isGeneric("printGenes"))
  setGeneric("printGenes", function(object, whichTerms, file, ...) standardGeneric("printGenes"))


########
## TODO
##
## remove the chip argument and replace it with a function which,
## given a list of genes will provide information about these genes
##
## TODO
########


## if the file argument is missing the function will just return
## a list of data.frames, each data.frame containg the gene information for specified GO term
setMethod("printGenes",
          signature(object = "topGOdata", whichTerms = "character", file = "missing"),
          function(object, whichTerms, chip, numChar = 100, simplify = TRUE,
                   geneCutOff = 50, pvalCutOff) {

            term.genes <- genesInTerm(object, whichTerms)
            all.genes <- character()
            lapply(term.genes, function(x) all.genes <<- union(x, all.genes))

            chip <- sub(".db$", "", chip)

            LL.lib <- get(paste(chip, "ENTREZID", sep = ""))
            Sym.lib <- get(paste(chip, "SYMBOL", sep = ""))
            GNAME.lib <- get(paste(chip, "GENENAME", sep = ""))

            probeMapping <- data.frame(LL.id = as.integer(unlist(mget(all.genes,
                                         envir = LL.lib, ifnotfound = NA))),
                                       Symbol.id = unlist(mget(all.genes,
                                         envir = Sym.lib, ifnotfound = NA)))

            retList <- vector("list", length(whichTerms))
            names(retList) <- whichTerms

            for(gt in whichTerms) {
              affID <- term.genes[[gt]]

              pval <- sort(geneScore(object, affID, use.names = TRUE))
              if(missing(pvalCutOff))
                affID <- names(pval)
              else {
                pval <- pval[pval <= pvalCutOff]
                affID <- names(pval)
              }

              ## we restrict the output to the number of genes
              length(affID) <- min(length(affID), geneCutOff)

              ## if there are no genes, there is nothing to print
              if(length(affID) == 0) {
                message("\n\t No genes over the cutOff for: ", gt)
                next
              }

              genesNames <- sapply(mget(affID, envir = GNAME.lib),
                                   function(x) return(x[1]))
              genesNames <- paste(substr(genesNames, 1, numChar),
                                  ifelse(nchar(genesNames) > numChar, "...", ""), sep = "")

              retList[[gt]] <- cbind("Chip ID" = affID, probeMapping[affID, ], "Gene name" = genesNames,
                                     "raw p-value" = format.pval(pval[affID], digits = 3, eps = 1e-30))
            }

            ## if we have only one GO term we return just the data.frame
            if(simplify && length(whichTerms) == 1)
              return(retList[[1]])

            return(retList)
          })


setMethod("printGenes",
          signature(object = "topGOdata", whichTerms = "character", file = "character"),
          ## numChar = "integer"),
          function(object, whichTerms, file, oneFile = FALSE, ...) {

            infoList <- printGenes(object, whichTerms, ...)

            if(length(whichTerms) == 1) {
              write.table(infoList, quote = TRUE, sep = ",", append = FALSE,
                          file = paste(paste(file, sub(":", "_", whichTerms), sep = "_"), "csv", sep = "."),
                          col.names = colnames(infoList), row.names = 1:nrow(infoList))
              return()
            }

            if(oneFile) {
              fAppend <- FALSE
              for(gt in whichTerms) {
                write.table(cbind(GOID = rep(gt, nrow(infoList[[gt]])), infoList[[gt]]),
                            quote = TRUE, sep = ",", append = fAppend,
                            file = paste(file, "csv", sep = "."),
                            col.names = colnames(infoList[[gt]]), row.names = FALSE)
                fAppend <- TRUE
                return()
              }
            }

            for(gt in whichTerms)
              write.table(infoList[[gt]], quote = TRUE, sep = ",", append = FALSE,
                          file = paste(paste(file, sub(":", "_", gt), sep = "_"), "csv", sep = "."),
                          col.names = colnames(infoList[[gt]]), row.names = 1:nrow(infoList[[gt]]))
          })

 ######################################################################


.getTermsDefinition <- function(whichTerms, ontology, numChar = 20, multipLines = FALSE) {

  qTerms <- paste(paste("'", whichTerms, "'", sep = ""), collapse = ",")
  retVal <- dbGetQuery(GO_dbconn(), paste("SELECT term, go_id FROM go_term WHERE ontology IN",
                                          "('", ontology, "') AND go_id IN (", qTerms, ");", sep = ""))

  termsNames <- retVal$term
  names(termsNames) <- retVal$go_id

  if(!multipLines)
    shortNames <- paste(substr(termsNames, 1, numChar),
                        ifelse(nchar(termsNames) > numChar, '...', ''), sep = '')
  else
    shortNames <- sapply(termsNames,
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\\\n"))
                         })

  names(shortNames) <- names(termsNames)
  
  ## return NAs for the terms that are not found in the DB and make sure the 'names' attribute is as specified
  shortNames <- shortNames[whichTerms]
  names(shortNames) <- whichTerms

  return(shortNames)
}


## methodsSig contains a named vector of p-values for each run method
.sigAllMethods <- function (methodsSig)
{
  names.index <- names(methodsSig[[1]])
  retval <- as.data.frame(lapply(methodsSig, function(x) x[names.index]))
  names(retval) <- names(methodsSig)
  return(retval)
}



######################################################################
if(!isGeneric("GenTable"))
  setGeneric("GenTable", function(object, ...) standardGeneric("GenTable"))

setMethod("GenTable",
          signature(object = "topGOdata"),
          ## ... = list of topGOresult object
          ## orderBy = "ANY", ## integer or character (index/name)
          ## ranksOf = "ANY", ## which ranks to be computed (integer/character)
          ## topNodes = "integer",
          ## numChar = "integer",
          ## useLevels = "logical"),
          function(object, ..., orderBy = 1, ranksOf = 2,
                   topNodes = 10, numChar = 40,
                   format.FUN = format.pval, decreasing = FALSE,
                   useLevels = FALSE) {

            resList <- list(...)

            ## first for the class of the elements in the list
            if(!all(sapply(resList, is, "topGOresult")))
              stop("Use: topGOdata, topGOresult_1, topGOresult_2, ..., \"parameters\".")

            ## if no names were provided we name them
            if(is.null(names(resList)))
              names(resList) <- paste("result", 1:length(resList), sep = "")

            ## obtain the score from the objects
            resList <- lapply(resList, score)

            ## order the scores and take care of the case in which only one result is provided
            ## in such case the orderBy and ranksOf parameters are ignored.
            if(length(resList) == 1) {
              orderBy <- ranksOf <- 1
              l <- data.frame(resList)
              names(l) <- ifelse(is.null(names(resList)), "", names(resList))
            } else {
              l <- .sigAllMethods(resList)
            }

            index <- order(l[, orderBy], decreasing = decreasing)
            l <- l[index, , drop = FALSE]

            if(decreasing)
              rr <- rank(-l[, ranksOf], ties.method = "first")
            else
              rr <- rank(l[, ranksOf], ties.method = "first")

            whichTerms <- rownames(l)[1:topNodes]
            l <- l[whichTerms, , drop = FALSE]
            rr <- as.integer(rr[1:topNodes])

            shortNames <- .getTermsDefinition(whichTerms, ontology(object), numChar = numChar)

            infoMat <- data.frame('GO ID' = whichTerms, 'Term' = shortNames, stringsAsFactors = FALSE)

            ## put the levels of the GO
            if(useLevels) {
              nodeLevel <- buildLevels(graph(object), leafs2root = TRUE)
              nodeLevel <- unlist(mget(whichTerms, envir = nodeLevel$nodes2level))
              infoMat <- data.frame(infoMat, Level = as.integer(nodeLevel))
            }

            annoStat <- termStat(object, whichTerms)

            ## if orderBy == ranksOf then there is no need to put the ranks
            if(ranksOf != orderBy) {
              dim(rr) <- c(length(rr), 1)
              colnames(rr) <- paste("Rank in ", ifelse(is.character(ranksOf), ranksOf, colnames(l)[ranksOf]), sep = "")

              infoMat <- data.frame(infoMat, annoStat, rr,
                                    apply(l, 2, format.FUN, dig = 2, eps = 1e-30),
                                    check.names = FALSE, stringsAsFactors = FALSE)

            } else {
              infoMat <- data.frame(infoMat, annoStat,
                                    apply(l, 2, format.FUN, dig = 2, eps = 1e-30),
                                    check.names = FALSE, stringsAsFactors = FALSE)
            }

            ##rownames(infoMat) <- whichTerms
            rownames(infoMat) <- 1:length(whichTerms)

            return(infoMat)
          })


## if(!isGeneric("genLatexTable"))
##   setGeneric("genLatexTable", function(object, resList, ...) standardGeneric("genLatexTable"))

## setMethod("genLatexTable",
##           signature(object = "topGOdata",
##                     resList = "list",
##                     orderBy = "ANY", ## integer or character (index/name)
##                     ranksOf = "ANY", ## which ranks to be computed (integer/character)
##                     topNodes = "integer",
##                     numChar = "integer"),
##           function(object, resList, orderBy, ranksOf, topNodes = 10, no.char = 40) {


## = infoMat,
##                         pval.tab = xtable.matrix(infoMat, caption = 'GO terms p-value', label = 'tab:GOinfo')))




################################################################################
## function that computes the raw/adjusted p-values based on
## multtest package

getPvalues <- function(edata, classlabel, test = "t",
                       alternative = c("greater", "two.sided", "less")[1],
                       genesID = NULL,
                       correction = c("none", "Bonferroni", "Holm", "Hochberg",
                         "SidakSS", "SidakSD", "BH", "BY")[8]) {

  require('multtest') || stop('package multtest is required')

  ## restrict the dataset
  if(!is.null(genesID))
    edata <- edata[intersect(genesID, rownames(edata)), ]
  genesID <- rownames(edata)

  t.stats <- mt.teststat(edata, classlabel = classlabel, test = test)

  if (alternative == "less")
    p.values <- pt(t.stats, df = length(classlabel) - 2)
  else if (alternative == "greater")
    p.values <- pt(t.stats, df = length(classlabel) - 2, lower.tail = FALSE)
  else
    p.values <- 2 * pt(-abs(t.stats), df = length(classlabel) - 2)

  if(correction != "none") {
    p.values <- mt.rawp2adjp(p.values, correction)
    p.values <- p.values$adjp[order(p.values$index), correction]
  }
  names(p.values) <- genesID

  return(p.values)
}

################################################################################

## a <= b
.sigRatio.ratio <- function(a, b, tolerance = 1e-50) {

  ## if a and b are almost equal we return 2
  if(identical(all.equal(a, b, tolerance = tolerance), TRUE))
    return(2)

  ## we want to compute b / a, thus we must take care for a = 0
  if(identical(all.equal(a, 0, tolerance = tolerance), TRUE))
    return(abs(b / (a + tolerance)))

  return(abs(b / a))
}

.sigRatio.log <- function(a, b, tolerance = 1e-50) {

  ## if a and b are almost equal we return 2
  if(identical(all.equal(a, b, tolerance = tolerance), TRUE))
    return(2)

  ## we want to compute log(a) / log(b), thus we must take care for b = 1
  if(identical(all.equal(log10(b), 0, tolerance = tolerance), TRUE))
    return(abs(log10(a) / tolerance))

  return(abs(log10(a) / log10(b)))
}

## a <= b
.sigRatio.01 <- function(a, b, tolerance = 1e-50) {

  ## if a and b are almost equal we return 2
  if(identical(all.equal(a, b, tolerance = tolerance), TRUE))
    return(2)

  if(a < b)
    return(1e50)

  return(2)
}

################################################################################
################################################################################

## functions to get gene stats from the topGOdata object
.getGeneData <- function(object) {
  return(c(Annotated = numGenes(object),
           Significant = numSigGenes(object),
           NodeSize = object@nodeSize))
}


## functions to get gene stats from the topGOdata object
.printGeneData <- function(x) {
  cat("Annotation data:\n")
  if("Annotated" %in% names(x))
    cat("    Annotated genes:", x["Annotated"], "\n")
  if("Significant" %in% names(x))
    cat("    Significant genes:", x["Significant"], "\n")
  cat("    Min. no. of genes annotated to a GO:", x["NodeSize"], "\n")
  if("SigTerms" %in% names(x))
    cat("    Nontrivial nodes:", x["SigTerms"], "\n")
}

################################################################################
################################################################################

## function to aggregate 2 or more topGOdata objects
combineResults <- function(..., method = c("gmean", "mean", "median", "min", "max")) {

  resList <- list(...)
  ## first for the class of the elements in the list
  if((length(resList) < 2) || !all(sapply(resList, is, "topGOresult")))
    stop("Use: topGOresult_1, topGOresult_2, ..., method = \"mean\"")

  combMethod <- match.arg(method)

  retVal <- resList[[1]]

  ## update the infos
  description(retVal) <- paste(description(retVal),
                               paste(combMethod, "(",
                                     paste(sapply(resList, algorithm), collapse = ", "),
                                     ")", sep = ""),
                               sep = "\n")
  testName(retVal) <- paste(unique(sapply(resList, testName)), collapse = ", ")
  algorithm(retVal) <- "ensemble"

  ## get the list of GOs
  nn <- names(score(retVal))
  ## obtain the score from the objects
  resList <- do.call(cbind, lapply(resList, score, whichGO = nn))

  newRes <- switch(combMethod,
                   mean = rowMeans(resList),
                   median = matrixStats::rowMedians(resList),
                   min = rowMins(resList),
                   max = rowMaxs(resList),
                   gmean = exp(rowMeans(log(resList))))
  names(newRes) <- rownames(resList) # just to make sure ... some of the functions we use might drop the names
  score(retVal) <- newRes

  return(retVal)
}
