#################### topGOdata methods ####################

##signature(.Object = "topGOdata",
##          ontology = "character",
##          allGenes = "ANY",
##          geneSelectionFun = "function",
##          description = "character"),
##          annotationFun = "function",

setMethod("initialize", "topGOdata",
          function(.Object,
                   ## which Ontology to be used
                   ontology,
                   ## a named numeric or factor, the names are the genes ID
                   allGenes,
                   ## function to select the signif. genes
                   geneSelectionFun = NULL,
                   ## description of the class
                   description = character(0),
                   ## expression matrix         #!
                   expressionMatrix = NULL,     #!
                   ## phenotype information     #!
                   phenotype = NULL,            #!
                   ## minimum node size
                   nodeSize = 1,
                   ## annotation function
                   annotationFun,
                   ## additional parameters for the annotationFun
                   ...) {

            .Object@description <- description
            .Object@ontology <- ontology

            .Object@expressionMatrix <- expressionMatrix  #!
            .Object@phenotype        <- phenotype         #!

            ## some checking
            if(is.null(names(allGenes)))
              stop("allGenes must be a named vector")

            if(!is.factor(allGenes) && !is.numeric(allGenes))
              stop("allGenes should be a factor or a numeric vector")

            .Object@allGenes <- names(allGenes)

            if(is.factor(allGenes)) {
              if(length(levels(allGenes)) != 2)
                stop("allGenes must be a factor with 2 levels")
              .Object@allScores <- factor(as.character(allGenes))
              .Object@geneSelectionFun <- function(x) {
                return(as.logical(as.integer(levels(x)))[x])
              }
            }
            else {
              .Object@allScores <- as.numeric(allGenes)

              ## function to select which genes are significant
              if(is.null(geneSelectionFun))
                warning("No function to select the significant genes provided!")
              .Object@geneSelectionFun <- geneSelectionFun
            }

            ## size of the nodes which will be pruned
            .Object@nodeSize = as.integer(max(nodeSize, 1))

            ## this function is returning a list of GO terms from the specified ontology
            ## whith each entry being a vector of genes
            message("\nBuilding most specific GOs .....")
            mostSpecificGOs <- annotationFun(ontology, .Object@allGenes, ...)
            message("\t( ", length(mostSpecificGOs), " GO terms found. )")

            ## the the GO graph is build started from the most specific terms
            message("\nBuild GO DAG topology ..........")
            g <- buildGOgraph.topology(names(mostSpecificGOs), ontology)
            message("\t( ",  numNodes(g), " GO terms and ", numEdges(g), " relations. )")

            ## probably is good to store the leves but for the moment we don't
            .nodeLevel <- buildLevels(g, leafs2root = TRUE)

            ## annotate the nodes in the GO graph with genes
            message("\nAnnotating nodes ...............")
            g <- mapGenes2GOgraph(g, mostSpecificGOs, nodeLevel = .nodeLevel) ## leafs2root

            ## select the feasible genes
            gRoot <- getGraphRoot(g)
            feasibleGenes <- ls(nodeData(g, n = gRoot, attr = "genes")[[gRoot]])
            message("\t( ", length(feasibleGenes), " genes annotated to the GO terms. )")

            .Object@feasible <- .Object@allGenes %in% feasibleGenes


            ## prune the GO graph
            if(.Object@nodeSize > 1) {
              cc <- .countsInNode(g, nodes(g))
              .Object@graph <- subGraph(names(cc)[cc >= .Object@nodeSize], g)
            } else {
              .Object@graph <-  g
            }

            .Object
          })



### temp  accessor methods ###
if(!isGeneric("expressionMatrix"))
  setGeneric("expressionMatrix", function(object) standardGeneric("expressionMatrix"))

setMethod("expressionMatrix", "topGOdata", function(object) object@expressionMatrix)

if(!isGeneric("phenotype"))
  setGeneric("phenotype", function(object) standardGeneric("phenotype"))

setMethod("phenotype", "topGOdata", function(object) object@phenotype)
######################################################################




#################### the accessor functions ####################

if(!isGeneric("description"))
  setGeneric("description", function(object) standardGeneric("description"))

setMethod("description", "topGOdata", function(object) object@description)

if(!isGeneric("ontology"))
  setGeneric("ontology", function(object) standardGeneric("ontology"))

setMethod("ontology", "topGOdata", function(object) object@ontology)

if(!isGeneric("allGenes"))
  setGeneric("allGenes", function(object) standardGeneric("allGenes"))

setMethod("allGenes", "topGOdata", function(object) object@allGenes)

##if(!isGeneric("allScores"))
##  setGeneric("allScores", function(object) standardGeneric("allScores"))

##setMethod("allScores", "topGOdata", function(object) object@allScores)

if(!isGeneric("feasible"))
  setGeneric("feasible", function(object) standardGeneric("feasible"))

setMethod("feasible", "topGOdata", function(object) object@feasible)

if(!isGeneric("graph"))
  setGeneric("graph", function(object) standardGeneric("graph"))

setMethod("graph", "topGOdata", function(object) object@graph)

## we take care of the nodes with fewer than "nodeSize" genes
##setMethod("graph", "topGOdata",
##           function(object) {
##             if(object@nodeSize > 0) {
##               cc <- .countsInNode(object@graph, nodes(object@graph))
##               return(subGraph(names(cc)[cc >= object@nodeSize], object@graph))
##             }
##             object@graph
##           })

if(!isGeneric("geneSelectionFun"))
  setGeneric("geneSelectionFun", function(object) standardGeneric("geneSelectionFun"))

setMethod("geneSelectionFun", "topGOdata", function(object) object@geneSelectionFun)


#################### the replacement functions ####################

if(!isGeneric("description<-"))
  setGeneric("description<-", function(object, value) standardGeneric("description<-"))

setMethod("description<-", "topGOdata", function(object, value) {object@description <- value; object})

if(!isGeneric("ontology<-"))
  setGeneric("ontology<-", function(object, value) standardGeneric("ontology<-"))

setMethod("ontology<-", "topGOdata", function(object, value) {object@ontology <- value; object})

## not a good idea to modify the list of allGenes, see updateGenes()
#if(!isGeneric("allGenes<-"))
#  setGeneric("allGenes<-", function(object, value) standardGeneric("allGenes<-"))
#
#setMethod("allGenes<-", "topGOdata", function(object, value) {object@allGenes <- value; object})


if(!isGeneric("feasible<-"))
  setGeneric("feasible<-", function(object, value) standardGeneric("feasible<-"))

setMethod("feasible<-", "topGOdata", function(object, value) {object@feasible <- value; object})


if(!isGeneric("geneSelectionFun<-"))
  setGeneric("geneSelectionFun<-", function(object, value) standardGeneric("geneSelectionFun<-"))

setMethod("geneSelectionFun<-", "topGOdata",
          function(object, value) {object@geneSelectionFun <- value; object})


if(!isGeneric("graph<-"))
  setGeneric("graph<-", function(object, value) standardGeneric("graph<-"))

setMethod("graph<-", "topGOdata", function(object, value) {object@graph <- value; object})



#################### methods to update/modify the objects ####################

## this function is used for updating the genes
## for the moment it allows to change the gene score and to restrict the gene names
## to the set of feasible genes
if(!isGeneric("updateGenes"))
  setGeneric("updateGenes", function(object, geneList, geneSelFun) standardGeneric("updateGenes"))

## for the case in which each gene has a score
setMethod("updateGenes",
          signature(object = "topGOdata", geneList = "numeric", geneSelFun = "function"),
          function(object, geneList, geneSelFun) {

            fGenes <- genes(object)

            if(is.null(names(geneList)))
              stop("geneList must be a named vector")

            if(!all(fGenes %in% names(geneList)))
              stop("Please provide a feasible geneList")

            object@allGenes <- names(geneList)

            object@allScores <- as.numeric(geneList)
            object@geneSelectionFun <- geneSelFun

            ## set up the feasible genes index vector
            object@feasible <- object@allGenes %in% fGenes

            ## TODO ........
            ## we need to update the graph
            ## object@graph <- updateGraph(object@graph, object@feasible)

            object
          })


setMethod("updateGenes",
          signature(object = "topGOdata", geneList = "factor", geneSelFun = "missing"),
          function(object, geneList) {

            fGenes <- genes(object)

            if(is.null(names(geneList)))
              stop("geneList must be a named vector")

            if(!all(genes(object) %in% names(geneList)))
              stop("Please provide a feasible geneList")

            object@allGenes <- names(geneList)

            if(length(levels(geneList)) != 2)
              stop("geneList must be a factor with 2 levels")
            object@allScores <- factor(as.character(geneList))
            object@geneSelectionFun <- function(x) {
              return(as.logical(as.integer(levels(x)))[x])
            }

            ## set up the feasible genes index vector
            object@feasible <- object@allGenes %in% fGenes

            ## TODO ........
            ## we need to update the graph
            ## object@graph <- updateGraph(object@graph, object@feasible)

            object
          })



#################### methods to obtain genes lists ####################

## return the genes that are used in the analysis (the one that can be
## mapped to the specific Ontology)
if(!isGeneric("genes"))
  setGeneric("genes", function(object) standardGeneric("genes"))

setMethod("genes", "topGOdata", function(object) object@allGenes[object@feasible])


if(!isGeneric("numGenes"))
  setGeneric("numGenes", function(object) standardGeneric("numGenes"))

setMethod("numGenes", "topGOdata", function(object) sum(object@feasible))


##
if(!isGeneric("geneScore"))
  setGeneric("geneScore", function(object, whichGenes, ...) standardGeneric("geneScore"))

setMethod("geneScore",
          signature(object = "topGOdata", whichGenes = "missing"),
          function(object, use.names = FALSE) {
            if(is.factor(object@allGenes))
              retList <- (as.numeric(levels(object@allScores))[object@allScores])[object@feasible]
            else
              retList <- as.numeric(object@allScores)[object@feasible]

            if(use.names)
              names(retList) <- genes(object)

            return(retList)
          })

setMethod("geneScore",
          signature(object = "topGOdata", whichGenes = "character"),
          function(object, whichGenes, use.names = TRUE) {
            index <- object@feasible & (object@allGenes %in% whichGenes)
            if(is.factor(object@allGenes))
              retList <- (as.numeric(levels(object@allScores))[object@allScores])[index]
            else
              retList <- as.numeric(object@allScores)[index]

            names(retList) <- object@allGenes[index]
            retList <- retList[whichGenes[whichGenes %in% genes(object)]]

            if(!use.names)
              names(retList) <- NULL

            return(retList)
          })


if(!isGeneric("sigGenes"))
  setGeneric("sigGenes", function(object) standardGeneric("sigGenes"))

setMethod("sigGenes", "topGOdata",
          function(object) {

            ##if(is.null(object@geneSelectionFun))
            ##  return(NULL)

            ## select the significant genes and the feasible ones
            sGenesIndex <- object@geneSelectionFun(object@allScores) & object@feasible
            return(object@allGenes[sGenesIndex])
          })

if(!isGeneric("numSigGenes"))
  setGeneric("numSigGenes", function(object) standardGeneric("numSigGenes"))

setMethod("numSigGenes", "topGOdata",
          function(object) {
            return(sum(object@geneSelectionFun(object@allScores) & object@feasible))
          })


############## methods to access the information on GO terms ##############
## we should exclusively use graph(object) to access the graph

if(!isGeneric("usedGO"))
  setGeneric("usedGO", function(object) standardGeneric("usedGO"))

setMethod("usedGO", "topGOdata", function(object) nodes(graph(object)))

if(!isGeneric("attrInTerm"))
  setGeneric("attrInTerm", function(object, attr, whichGO) standardGeneric("attrInTerm"))

setMethod("attrInTerm",
          signature(object = "topGOdata", attr = "character", whichGO = "character"),
          function(object, attr, whichGO) {
            return(.getFromNode(graph(object), attr, whichGO))
          })

setMethod("attrInTerm",
          signature(object = "topGOdata", attr = "character", whichGO = "missing"),
          function(object, attr, whichGO) {
            return(.getFromNode(graph(object), attr, nodes(graph(object))))
          })

## function that return for each GO term specified in the whichGO parameter the vector
## of annotated genes. If the whichGO is missing than the vector of annotated genes is
## returned for each GO term
if(!isGeneric("genesInTerm"))
  setGeneric("genesInTerm", function(object, whichGO) standardGeneric("genesInTerm"))

setMethod("genesInTerm",
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO) {
            return(.genesInNode(graph(object), whichGO))
          })


setMethod("genesInTerm",
          signature(object = "topGOdata", whichGO = "missing"),
          function(object) {
            return(.genesInNode(graph(object), nodes(graph(object))))
          })

## similar like above but returns the scores in each GO
if(!isGeneric("scoresInTerm"))
  setGeneric("scoresInTerm", function(object, whichGO, ...) standardGeneric("scoresInTerm"))

setMethod("scoresInTerm",
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO, use.names = FALSE) {
            l <- lapply(.genesInNode(graph(object), whichGO),
                        function(x) geneScore(object, x, use.names = use.names))
            return(l)
          })


setMethod("scoresInTerm",
          signature(object = "topGOdata", whichGO = "missing"),
          function(object, use.names = FALSE) {
            return(scoreInNode(object, nodes(graph(object)), use.names = use.names))
          })



if(!isGeneric("countGenesInTerm"))
  setGeneric("countGenesInTerm", function(object, whichGO) standardGeneric("countGenesInTerm"))

setMethod("countGenesInTerm",
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO) {
            return(.countsInNode(graph(object), whichGO))
          })

setMethod("countGenesInTerm",
          signature(object = "topGOdata", whichGO = "missing"),
          function(object) {
            return(.countsInNode(graph(object), nodes(graph(object))))
          })

## function that return for each GO term specified in the whichGO parameter
## a dataframe withe the following informations:
##   annotated - the number of annotated genes
##   significant  - the number of annotated significant genes
##   expected - how many genes are expected in a random case
if(!isGeneric("termStat"))
  setGeneric("termStat", function(object, whichGO) standardGeneric("termStat"))

setMethod("termStat",
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO) {

			if(length(whichGO) == 0) {
				return(data.frame(Annotated = integer(),
                             	  Significant = integer(),
                              	  Expected = numeric()))
			}

            x <- .genesInNode(graph(object), whichGO)

            anno <- sapply(x, length)

            sGenes <- sigGenes(object)
            sig <- sapply(x, function(e) length(intersect(e, sGenes)))

            expect <- numSigGenes(object) / length(genes(object))
            expect <- round(anno * expect, 2)

            return(data.frame(Annotated = anno,
                              Significant = sig,
                              Expected = expect,
                              row.names = names(x)))
          })

setMethod("termStat",
          signature(object = "topGOdata", whichGO = "missing"),
          function(object) termStat(object, nodes(graph(object))))



## write information in the graph nodes
if(!isGeneric("updateTerm<-"))
  setGeneric("updateTerm<-", function(object, attr, value)  standardGeneric("updateTerm<-"))

setMethod("updateTerm<-",
          signature(object = "topGOdata", attr = "character", value = "ANY"),
          function(object, attr, value) {
            ## we need to update the full graph ....
            ## object@graph <- .writeToNodes(graph(object), attr, value)
            object@graph <- .writeToNodes(object@graph, attr, value)
            object
          })


##############################  printing   ##############################
## make use of both "print" and "show" generics


setMethod("print", "topGOdata", function(x, ...) .printTopGOdata(x))
setMethod("show", "topGOdata", function(object) .printTopGOdata(x = object))

.printTopGOdata <- function(x) {
  cat("\n------------------------- topGOdata object -------------------------\n")
  cat("\n Description:\n")
  cat("   - ", x@description, "\n")

  cat("\n Ontology:\n")
  cat("   - ", x@ontology, "\n")

  ## all genes from the array
  cat("\n", length(x@allGenes) ,"available genes (all genes from the array):\n")
  sym <- x@allGenes[1:min(length(x@allGenes), 5)]
  cat("   - symbol: ", sym, " ...\n")
  if(is.numeric(x@allScores)) {
    score <- x@allScores[1:min(length(x@allGenes), 5)]
    score <- apply(cbind(score, nchar(sym)), 1, function(x) format(x[1], digits = max(x[2] - 2, 1)))
    cat("   - score : ", score, " ...\n")
  }
  cat("   -", sum(x@geneSelectionFun(x@allScores)), " significant genes. \n")

  ## feasible genes
  cat("\n", numGenes(x) ,"feasible genes (genes that can be used in the analysis):\n")
  sym <- genes(x)[1:min(numGenes(x), 5)]
  cat("   - symbol: ", sym, " ...\n")
  if(is.numeric(x@allScores)) {
    score <- geneScore(x)[1:min(numGenes(x), 5)]
    score <- apply(cbind(score, nchar(sym)), 1, function(x) format(x[1], digits = max(x[2] - 2, 1)))
    cat("   - score : ", score, " ...\n")
  }
  cat("   -", numSigGenes(x), " significant genes. \n")

  cat("\n GO graph (nodes with at least ", max(x@nodeSize, 1), " genes):\n")
  gg <- graph(x)
  cat("   - a graph with", edgemode(gg), "edges\n")
  cat("   - number of nodes =", numNodes(gg), "\n")
  cat("   - number of edges =", numEdges(gg), "\n")

  ##cat("\n Signif. genes sellection:\n\n")
  ##print(x@geneSelectionFun)
  cat("\n------------------------- topGOdata object -------------------------\n\n")
}



## TODO ....
## function to return the name of GO terms which have the no. of annotated
## genes in some interval


######################################################################
######################################################################


######################## topGOresult methods ########################

setMethod("initialize", "topGOresult",
          function(.Object, description = character(),
                   score, testName, algorithm, geneData = integer()) {
            .Object@description <- description
            .Object@score <- score
            .Object@testName <- testName
            .Object@algorithm <- algorithm
            .Object@geneData <- geneData

            .Object
          })


#################### the accessor functions ####################

setMethod("description", "topGOresult", function(object) object@description)

setMethod("score", "topGOresult",
          function(x, whichGO, ...) {

            if(missing(whichGO))
              return(x@score)

            if(!is.character(whichGO))
              stop("whichGO must be a valid GO id!")

            allGO <- names(x@score)
            feasableGO <- intersect(whichGO, allGO) # in this way the order of whichGO is preserved

            if(length(feasableGO) < 1L) {
              warning("The specified GO terms could not be found in the object... Returning empty vector")
              return(numeric())
            }

            if(length(setdiff(whichGO, allGO)) > 0L)
              warning("Not all specified GOs can be found in the object.")
            if(length(setdiff(allGO, whichGO)) > 0L)
              warning("Not all GOs from the object were retrived.")

            return(x@score[feasableGO])
          })

if(!isGeneric("testName"))
  setGeneric("testName", function(object) standardGeneric("testName"))

setMethod("testName", "topGOresult", function(object) object@testName)


if(!isGeneric("algorithm"))
  setGeneric("algorithm", function(object) standardGeneric("algorithm"))

setMethod("algorithm", "topGOresult", function(object) object@algorithm)


if(!isGeneric("geneData"))
  setGeneric("geneData", function(object) standardGeneric("geneData"))

setMethod("geneData", "topGOresult", function(object) object@geneData)

#################### the replacement functions ####################

setMethod("description<-", "topGOresult",
          function(object, value) {object@description <- value; object})

setMethod("score<-", "topGOresult",
          function(x, ..., value) {
            x@score <- value
            x
          })


if(!isGeneric("testName<-"))
  setGeneric("testName<-", function(object, value) standardGeneric("testName<-"))

setMethod("testName<-", "topGOresult", function(object, value) {object@testName <- value; object})


if(!isGeneric("algorithm<-"))
  setGeneric("algorithm<-", function(object, value) standardGeneric("algorithm<-"))

setMethod("algorithm<-", "topGOresult", function(object, value) {object@algorithm <- value; object})


if(!isGeneric("geneData<-"))
  setGeneric("geneData<-", function(object, value) standardGeneric("geneData<-"))

setMethod("geneData<-", "topGOresult", function(object, value) {object@geneData <- value; object})


##############################  printing   ##############################
## make use of both "print" and "show" generics

setMethod("print", "topGOresult", function(x, ...) .printTopGOresult(x))
setMethod("show", "topGOresult", function(object) .printTopGOresult(x = object))

.printTopGOresult <- function(x) {
  cat("\nDescription:", description(x), "\n")
  cat("'", algorithm(x), "' algorithm with the '", testName(x), "' test\n", sep = "")
  cat(length(score(x)), "GO terms scored:", sum(score(x) <=  0.01), "terms with p < 0.01\n")
  .printGeneData(geneData(x))
}





######################################################################
######################################################################



######################## groupStats methods ########################

setMethod("initialize", "groupStats",
          function(.Object, testStatistic, name, allMembers, groupMembers) {
            .Object@name <- name
            .Object@allMembers <- allMembers
            .Object@members <- groupMembers
            .Object@testStatistic <- testStatistic
            ##.Object@testStatPar <- testStatPar

            .Object
          })


#################### the accessor functions ####################

if(!isGeneric("Name"))
  setGeneric("Name", function(object) standardGeneric("Name"))

setMethod("Name", "groupStats", function(object) object@name)

if(!isGeneric("allMembers"))
  setGeneric("allMembers", function(object) standardGeneric("allMembers"))

setMethod("allMembers", "groupStats", function(object) object@allMembers)

##if(!isGeneric("members"))
setGeneric("members", function(x, i) standardGeneric("members"))

setMethod("members",
          signature(x = "groupStats", i = "missing"),
          function(x, i) x@members)

if(!isGeneric("testStatistic"))
  setGeneric("testStatistic", function(object) standardGeneric("testStatistic"))

setMethod("testStatistic", "groupStats", function(object) object@testStatistic)

if(!isGeneric("testStatPar"))
  setGeneric("testStatPar", function(object) standardGeneric("testStatPar"))

setMethod("testStatPar", "groupStats", function(object) object@testStatPar)

#################### the replacement functions ####################

if(!isGeneric("Name<-"))
  setGeneric("Name<-", function(object, value) standardGeneric("Name<-"))

setMethod("Name<-", "groupStats", function(object, value) {object@Name <- value; object})

if(!isGeneric("allMembers<-"))
  setGeneric("allMembers<-", function(object, value) standardGeneric("allMembers<-"))

setMethod("allMembers<-", "groupStats", function(object, value) {object@allMembers <- value; object})

if(!isGeneric("members<-"))
  setGeneric("members<-", function(object, value) standardGeneric("members<-"))

setMethod("members<-", "groupStats", function(object, value) {object@members <- value; object})


#################### other functions ####################

if(!isGeneric("numMembers"))
  setGeneric("numMembers", function(object) standardGeneric("numMembers"))

setMethod("numMembers", "groupStats", function(object) length(object@members))

if(!isGeneric("numAllMembers"))
  setGeneric("numAllMembers", function(object) standardGeneric("numAllMembers"))

setMethod("numAllMembers", "groupStats", function(object) length(object@allMembers))


## MAIN function -- it should return the "p-value" of the Test Satatistic
if(!isGeneric("runTest"))
  setGeneric("runTest", function(object, algorithm, statistic, ...) standardGeneric("runTest"))

setMethod("runTest",
          signature(object = "groupStats", algorithm = "missing", statistic = "missing"),
          function(object) object@testStatistic(object))


## function to update/(build a new object)
if(!isGeneric("updateGroup"))
  setGeneric("updateGroup", function(object, name, members, ...) standardGeneric("updateGroup"))

setMethod("updateGroup",
          signature(object = "groupStats", name = "character", members = "character"),
          function(object, name, members) {
            object@name <- name
            object@members <- members

            return(object)
          })



######################################################################
######################################################################


#################### classicCount methods ####################

#################### constructor ####################
setMethod("initialize", "classicCount",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers)
            .Object@significant <- which(allMembers %in% sigMembers)
            .Object@testStatPar = list(...)

            .Object
          })


if(!isGeneric("sigMembers<-"))
  setGeneric("sigMembers<-", function(object, value) standardGeneric("sigMembers<-"))

setMethod("sigMembers<-", "classicCount",
          function(object, value) {
            object@significant <- which(object@allMembers %in% value)
            object
          })

if(!isGeneric("sigAllMembers"))
  setGeneric("sigAllMembers", function(object) standardGeneric("sigAllMembers"))

setMethod("sigAllMembers", "classicCount",
          function(object) object@allMembers[object@significant])

if(!isGeneric("numSigAll"))
  setGeneric("numSigAll", function(object) standardGeneric("numSigAll"))

setMethod("numSigAll", "classicCount",
          function(object) length(object@significant))

if(!isGeneric("sigMembers"))
  setGeneric("sigMembers", function(object) standardGeneric("sigMembers"))

setMethod("sigMembers", "classicCount",
          function(object) intersect(sigAllMembers(object), members(object)))

if(!isGeneric("numSigMembers"))
  setGeneric("numSigMembers", function(object) standardGeneric("numSigMembers"))

setMethod("numSigMembers", "classicCount",
          function(object) sum(members(object) %in% sigAllMembers(object)))

if(!isGeneric("contTable"))
  setGeneric("contTable", function(object) standardGeneric("contTable"))

setMethod("contTable", "classicCount",
          function(object) {

            numHits <- numSigMembers(object)
            numSig <- numSigAll(object)
            numM <- numMembers(object)

            contMat <- cbind(sig = c(numHits, numSig - numHits),
                             notSig = c(numM - numHits, numAllMembers(object) - numM - numSig + numHits))
            row.names(contMat) <- c("anno", "notAnno")

            return(contMat)
          })



#################### classicScore methods ####################

setMethod("initialize", "classicScore",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   score = numeric(),
                   scoreOrder = "increasing",  ## this is the case in which p-values are used
                   ...) {

            scoreOrder <- switch(scoreOrder,
                                 decreasing = TRUE,
                                 increasing = FALSE,
                                 stop("scoreOrder should be increasing or decreasing"))

            ## first we order the members according to the score
            index <- order(score, decreasing = scoreOrder)
            if(length(allMembers) != length(score))
              warning("score length don't match.")

            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers[index], groupMembers)

            .Object@score <- as.numeric(score)[index]
            .Object@scoreOrder <- scoreOrder
            .Object@testStatPar = list(...)

            .Object
          })

if(!isGeneric("scoreOrder"))
  setGeneric("scoreOrder", function(object) standardGeneric("scoreOrder"))

setMethod("scoreOrder", "classicScore", function(object) object@scoreOrder)

## methods to get the score
if(!isGeneric("allScore"))
  setGeneric("allScore", function(object, use.names) standardGeneric("allScore"))

setMethod("allScore",
          signature(object = "classicScore", use.names = "missing"),
          function(object) object@score)

setMethod("allScore",
          signature(object = "classicScore", use.names = "logical"),
          function(object, use.names = FALSE) {
            if(use.names) {
              ret.val <- object@score
              names(ret.val) <- allMembers(object)
              return(ret.val)
            }

            object@score
          })


if(!isGeneric("membersScore"))
  setGeneric("membersScore", function(object) standardGeneric("membersScore"))

setMethod("membersScore", "classicScore",
          function(object) {
            index <- allMembers(object) %in% members(object)
            ss <- object@score[index]
            names(ss) <- allMembers(object)[index]

            ## just to be safe of the order
            return(ss[members(object)])
          })


## methods to assign the score
## the value should be a named vector, the names should be allMembers
setMethod("score<-", "classicScore",
          function(x, ..., value) {
            ## some checking
            if(length(names(value)) == 0)
              stop("no names associated with the score")
            if(!all(names(value) %in% x@allMembers))
              warning("The new score names do not match with allMembers.")
            if(!all(x@members %in% names(value)))
              stop("You need to build a new object!")

            value <- sort(value, x@scoreOrder)
            x@allMembers <- names(value)
            x@score <- as.numeric(value)

            return(x)
          })

## rank the members of the group according to their score
if(!isGeneric("rankMembers"))
  setGeneric("rankMembers", function(object) standardGeneric("rankMembers"))

## return(which(object@allMembers %in% object@members)) replace with match ... should be faster
setMethod("rankMembers", "classicScore",
          function(object) {
            return(match(object@members, object@allMembers))
          })



#################### classicExpr methods ####################

#################### constructor ####################
setMethod("initialize", "classicExpr",
          function(.Object,
                   testStatistic,
                   name = character(),
                   ## allMembers = character(), ## this should be given in the eData matrix
                   groupMembers = character(),
                   exprDat,
                   pType = factor(),
                   ...) {

            if(missing(exprDat)) {
              allMembers <- character()
              e <- emptyenv()
            }
            else {
              if(class(exprDat) != "matrix")
                error("exprDat must be of type matrix")

              allMembers <- rownames(exprDat)
              e <- new.env(hash = TRUE, parent = emptyenv())
              rownames(exprDat) <- NULL
              assign("expresMatrix", exprDat, envir = e)
            }

            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers)

            .Object@pType <- pType
            .Object@eData <- e
            .Object@testStatPar = list(...)

            .Object
          })


if(!isGeneric("pType<-"))
  setGeneric("pType<-", function(object, value) standardGeneric("pType<-"))

setMethod("pType<-", "classicExpr", function(object, value) {object@pType <- value; object})

if(!isGeneric("pType"))
  setGeneric("pType", function(object) standardGeneric("pType"))

## this is for the simple case in which the pType is a vector!
setMethod("pType", "classicExpr", function(object) object@pType)

setMethod("allMembers<-", "classicExpr",
           function(object, value) {warning("Assignmanet not allowed"); object})

if(!isGeneric("emptyExpr"))
  setGeneric("emptyExpr", function(object) standardGeneric("emptyExpr"))

setMethod("emptyExpr", "classicExpr",
          function(object) return(environmentName(object@eData) == "R_EmptyEnv"))

if(!isGeneric("membersExpr"))
  setGeneric("membersExpr", function(object) standardGeneric("membersExpr"))

setMethod("membersExpr", "classicExpr", function(object)
          get("expresMatrix", envir = object@eData)[object@allMembers %in% members(object), , drop = FALSE])



######################################################################
######################################################################


#################### weight01Count methods ####################

#################### constructor ####################
setMethod("initialize", "weight01Count",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   elim = integer(),
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers,
                                      sigMembers)

            .Object@elim <- which(.Object@members %in% elim)
            .Object@testStatPar = list(...)

            .Object
          })

if(!isGeneric("elim<-"))
  setGeneric("elim<-", function(object, value) standardGeneric("elim<-"))

setMethod("elim<-", "weight01Count",
          function(object, value) {
            object@elim <- which(object@members %in% value)
            object
          })

if(!isGeneric("elim"))
  setGeneric("elim", function(object) standardGeneric("elim"))

setMethod("elim", "weight01Count",
          function(object) object@members[object@elim])

## account for the eliminated members
#TODO:

setMethod("numMembers", "weight01Count",
          function(object) length(object@members) - length(object@elim))

setMethod("numAllMembers", "weight01Count",
          function(object) length(object@allMembers) - length(object@elim))

setMethod("sigAllMembers", "weight01Count",
          function(object) setdiff(object@allMembers[object@significant], elim(object)))

setMethod("numSigAll", "weight01Count",
          function(object) length(sigAllMembers(object)))

setMethod("sigMembers", "weight01Count",
          function(object) intersect(sigAllMembers(object), members(object)))

setMethod("numSigMembers", "weight01Count",
          function(object) length(sigMembers(object)))





#################### weight01Score methods ####################

#################### constructor ####################
setMethod("initialize", "weight01Score",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   score = numeric(),
                   scoreOrder = "increasing",
                   elim = integer(),
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers, score,
                                      scoreOrder)

            .Object@elim <- which(.Object@members %in% elim)
            .Object@testStatPar = list(...)

            .Object
          })

setMethod("elim<-", "weight01Score",
          function(object, value) {
            object@elim <- which(object@members %in% value)
            object
          })

setMethod("elim", "weight01Score", function(object) object@members[object@elim])

setMethod("members",
          signature(x = "weight01Score", i = "missing"),
          function(x) {
            if(length(x@elim) == 0)
              return(x@members)

            return(x@members[-x@elim])
          })

setMethod("allMembers", "weight01Score",
          function(object) object@allMembers[!(object@allMembers %in% elim(object))])

setMethod("numMembers", "weight01Score",
          function(object) length(object@members) - length(object@elim))

setMethod("numAllMembers", "weight01Score",
          function(object) length(object@allMembers) - length(object@elim))

setMethod("allScore",
          signature(object = "weight01Score", use.names = "missing"),
          function(object) {
            ret.val <- object@score
            index <- !(object@allMembers %in% elim(object))

            return(ret.val[index])
          })

setMethod("allScore",
          signature(object = "weight01Score", use.names = "logical"),
          function(object, use.names = FALSE) {
            ret.val <- object@score
            index <- !(object@allMembers %in% elim(object))

            if(use.names)
              names(ret.val) <- object@allMembers

            return(ret.val[index])
          })

setMethod("membersScore", "weight01Score",
          function(object) {
            ## since members(object) has the eliminated members removed
            index <- object@allMembers %in% members(object)

            ss <- object@score[index]
            names(ss) <- object@allMembers[index]

            ## just to be safe of the order
            return(ss[members(object)])
          })


setMethod("rankMembers", "weight01Score",
          function(object) {
            ## NEED to use methods allMembers and members since they accont for removed members
            return(which(allMembers(object) %in% members(object)))
          })




#################### weight01Expr methods ####################

#################### constructor ####################
setMethod("initialize", "weight01Expr",
          function(.Object,
                   testStatistic,
                   name = character(),
                   groupMembers = character(),
                   exprDat,
                   pType = factor(),
                   elim = integer(),
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      groupMembers, exprDat, pType)

            .Object@elim <- which(.Object@members %in% elim)
            .Object@testStatPar = list(...)

            .Object
          })

if(!isGeneric("elim<-"))
  setGeneric("elim<-", function(object, value) standardGeneric("elim<-"))

setMethod("elim<-", "weight01Expr",
          function(object, value) {
            object@elim <- which(object@members %in% value)
            object
          })

if(!isGeneric("elim"))
  setGeneric("elim", function(object) standardGeneric("elim"))

setMethod("elim", "weight01Expr",
          function(object) object@members[object@elim])

## probably not a good idea ....
setMethod("members",
          signature(x = "weight01Expr", i = "missing"),
          function(x) {
            if(length(x@elim) == 0)
              return(x@members)

            return(x@members[-x@elim])
          })

setMethod("allMembers", "weight01Expr",
          function(object) callNextMethod(object)[!(callNextMethod(object) %in% elim(object))])

setMethod("numMembers", "weight01Expr",
          function(object) length(object@members) - length(object@elim))

setMethod("numAllMembers", "weight01Expr",
          function(object) length(object@allMembers) - length(object@elim))



######################################################################
######################################################################




#################### elimCount methods ####################

#################### constructor ####################
setMethod("initialize", "elimCount",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   elim = integer(),
                   cutOff = 0.01,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers,
                                      sigMembers, elim = elim)
            .Object@cutOff <- cutOff
            .Object@testStatPar = list(...)

            .Object
          })

if(!isGeneric("cutOff<-"))
  setGeneric("cutOff<-", function(object, value) standardGeneric("cutOff<-"))

setMethod("cutOff<-", "elimCount",
          function(object, value)  {object@cutOff <- value; object})


if(!isGeneric("cutOff"))
  setGeneric("cutOff", function(object) standardGeneric("cutOff"))

setMethod("cutOff", "elimCount", function(object) object@cutOff)



#################### elimScore methods ####################

#################### constructor ####################
setMethod("initialize", "elimScore",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   score = numeric(),
                   scoreOrder = "increasing",
                   elim = integer(),
                   cutOff = 0.01,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers, score,
                                      scoreOrder, elim = elim)
            .Object@cutOff <- cutOff
            .Object@testStatPar = list(...)

            .Object
          })

setMethod("cutOff<-", "elimScore",
          function(object, value)  {object@cutOff <- value; object})

setMethod("cutOff", "elimScore", function(object) object@cutOff)



#################### elimExpr methods ####################

#################### constructor ####################
setMethod("initialize", "elimExpr",
          function(.Object,
                   testStatistic,
                   name = character(),
                   groupMembers = character(),
                   exprDat,
                   pType = factor(),
                   elim = integer(),
                   cutOff = 0.01,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      groupMembers, exprDat, pType,
                                      elim = elim)
            .Object@cutOff <- cutOff
            .Object@testStatPar = list(...)

            .Object
          })

setMethod("cutOff<-", "elimExpr",
          function(object, value)  {object@cutOff <- value; object})

setMethod("cutOff", "elimExpr", function(object) object@cutOff)



######################################################################
######################################################################



#################### weightCount methods ####################

#################### constructor ####################
setMethod("initialize", "weightCount",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   weights = numeric(),
                   sigRatio = "ratio",
                   penalise = "more",
                   ...) {

            .Object <- callNextMethod(.Object, testStatistic,
                                      paste(name, sigRatio, sep = " : "),
                                      allMembers, groupMembers, sigMembers)

            if(length(weights) == 0)
              .Object@weights <- numeric()
            else
              if(is.null(names(weights)) &&
                 length(intersect(names(weights), .Object@members)) != length(.Object@members))
                stop("The weight vector must be a named vector.")

            .Object@weights <- as.numeric(weights[.Object@members]) ## to remove the names

            ## probably the easiest way to store which sigRatio function is used
            .Object@sigRatio <- switch(sigRatio,
                                       ratio = .sigRatio.ratio,
                                       log = .sigRatio.log,
                                       "01" = .sigRatio.01,
                                       stop("Please give a valid sigRation name"))

            if(penalise == "more")
              ## need to think about this more
              .Object@penalise <- function(a, b) return(1 / (abs(log(a * b) / 2) + 1))
            else
              .Object@penalise <- function(a, b) return(1)

            .Object@roundFun <- floor
            .Object@testStatPar = list(...)

            .Object
          })

if(!isGeneric("penalise"))
  setGeneric("penalise", function(object, a, b) standardGeneric("penalise"))

setMethod("penalise",
          signature(object = "weightCount", a = "numeric", b = "numeric"),
          function(object, a, b) object@penalise(a, b))


setMethod("updateGroup",
          signature(object = "weightCount", name = "character", members = "character"),
          function(object, name, members, weights) {
            object <- callNextMethod(object, name, members)
            if(!missing(weights) && length(members) > 0) {
              if(length(intersect(names(weights), members)) != length(members))
                stop("weights and members vectors don't agree.")

              object@weights <- as.numeric(weights[members])
            }

            return(object)
          })

if(!isGeneric("Weights"))
  setGeneric("Weights", function(object, use.names) standardGeneric("Weights"))

setMethod("Weights",
          signature(object = "weightCount", use.names = "missing"),
          function(object) object@weights)

setMethod("Weights",
          signature(object = "weightCount", use.names = "logical"),
          function(object, use.names = FALSE) {
            if(use.names) {
              ret.val <- object@weights
              names(ret.val) <- members(object)
              return(ret.val)
            }

            object@weights
          })


if(!isGeneric("Weights<-"))
  setGeneric("Weights<-", function(object, value) standardGeneric("Weights<-"))

setMethod("Weights<-", "weightCount",
          function(object, value) {
            object@weights <- as.numeric(value[object@members])
            object
          })


#if(!isGeneric("sigRatio"))
#  setGeneric("sigRatio", function(object) standardGeneric("sigRatio"))
#
#setMethod("sigRatio", "weightCount", function(object) object@sigRatio)

if(!isGeneric("sigRatio<-"))
  setGeneric("sigRatio<-", function(object, value) standardGeneric("sigRatio<-"))

setMethod("sigRatio<-", "weightCount",
          function(object, value)  {
            object@sigRatio <- switch(value,
                                      log = .sigRatio.log,
                                      ratio = .sigRatio.ratio,
                                      "01" = .sigRatio.01,
                                      stop("Please give a valid sigRation name"))
            object
          })


if(!isGeneric("getSigRatio"))
  setGeneric("getSigRatio", function(object, a, b) standardGeneric("getSigRatio"))

setMethod("getSigRatio", "weightCount",
          function(object, a, b) object@sigRatio(a, b))


## account for the weights of the members
#TODO:

setMethod("numMembers", "weightCount",
          function(object) object@roundFun(sum(object@weights)))

setMethod("numAllMembers", "weightCount",
          function(object) {
            a <- sum(object@weights) + (length(object@allMembers) - length(object@members))
            return(object@roundFun(a))
          })

setMethod("numSigAll", "weightCount",
          function(object) {
            ## which group members are sig.
            index <- members(object) %in% sigAllMembers(object)
            ## how many sig. members in total
            num.ones <- length(object@significant) - sum(index)

            return(object@roundFun(sum(object@weights[index]) + num.ones))
          })

setMethod("numSigMembers", "weightCount",
          function(object) {
            ## sum the weights of the group members which are sig.
            return(object@roundFun(sum(object@weights[object@members %in% sigAllMembers(object)])))
          })



######################################################################
######################################################################

#################### parentChild methods ####################

#################### constructor ####################
setMethod("initialize", "parentChild",
          function(.Object,
                   testStatistic,
                   name = character(),
                   groupMembers = character(),
                   parents, ## should be a list
                   sigMembers = character(),
                   joinFun = c("union", "intersect"),
                   ...) {

            if(missing(parents)) {
              allMembers <- character()
              splitIndex <- integer()
              ## no use to have significant members if we don't
              ## know all the members
              sigMembers <- integer()
            }
            else {
              splitIndex <- sapply(parents, length)
              allMembers <- unlist(parents, use.names = FALSE)
              ## if some elements of sigMembers are not in allMembers
              sigMembers <- match(sigMembers, allMembers)
              sigMembers <- sigMembers[!is.na(sigMembers)]
            }

            .Object <- callNextMethod(.Object, testStatistic = testStatistic,
                                      name = name, allMembers = allMembers,
                                      groupMembers = groupMembers)

            .Object@joinFun <- match.arg(joinFun)
            .Object@significant <- sigMembers
            .Object@splitIndex <- splitIndex
            .Object@testStatPar = list(...)

            .Object
          })


if(!isGeneric("joinFun"))
  setGeneric("joinFun", function(object) standardGeneric("joinFun"))

setMethod("joinFun", "parentChild", function(object) object@joinFun)


## the main function is "allMembers"
## it selects the "parent" members depending on the joinFun function
setMethod("allMembers", "parentChild",
          function(object) {
            if(joinFun(object) == "union")
              return(unique(object@allMembers))

            ## split the vector into a list
            ss <- rep(1:length(object@splitIndex), times = object@splitIndex)
            l <- split(object@allMembers, ss)

            ## works for a list with at least one element
            res <- l[[1]]
            for(s in l[-1]) res <- intersect(s, res)

            return(res)
          })


if(!isGeneric("allParents"))
  setGeneric("allParents", function(object) standardGeneric("allParents"))

setMethod("allParents", "parentChild",
          function(object) {
            ## split the vector into a list
            ss <- rep(1:length(object@splitIndex), times = object@splitIndex)
            l <- object@allMembers
            names(l) <- names(ss)
            return(l)
          })

## needs to be redifined such that it accounts for the joinFun
setMethod("numAllMembers", "parentChild", function(object) length(allMembers(object)))

setMethod("sigAllMembers", "parentChild",
          function(object) intersect(object@allMembers[object@significant], allMembers(object)))

setMethod("numSigAll", "parentChild",
          function(object) length(sigAllMembers(object)))


setMethod("sigMembers<-", "parentChild",
          function(object, value) {
            sig <- match(value, object@allMembers)
            object@significant <- sig[!is.na(sig)]
            object
          })

setMethod("allMembers<-", "parentChild",
          function(object, value) {warning("Assignmanet not allowed"); object})

setMethod("updateGroup",
          signature(object = "parentChild", name = "missing", members = "character"),
          function(object, members, parents, sigMembers) {

            object@splitIndex <- sapply(parents, length)
            object@allMembers <- unlist(parents, use.names = FALSE)

            sigMembers <- match(sigMembers, object@allMembers)
            object@significant <- sigMembers[!is.na(sigMembers)]

            object@members <- members

            return(object)
          })





#################### pC methods ####################

#################### constructor ####################
setMethod("initialize", "pC",
          function(.Object,
                   testStatistic,
                   name = character(),
                   groupMembers = character(),
                   parents, ## should be a list
                   sigMembers = character(),
                   joinFun = c("union", "intersect"),
                   ...) {

            joinFun <- match.arg(joinFun)
            name <- paste(name, paste("joinFun = ", joinFun, sep = ""), sep = " : ")

            joinFun <- switch(joinFun,
                              "union" = get("union"),
                              "intersect" = get("intersect"))

            if(missing(parents)) {
              allMembers <- character()
              sigMembers <- integer()
            }
            else {
              ## works for a list with at least one element
              allMembers <- parents[[1]]
              for(s in parents[-1]) allMembers <- joinFun(s, allMembers)
              ## if some elements of sigMembers are not in allMembers
              sigMembers <- which(allMembers %in% sigMembers)
            }

            .Object <- callNextMethod(.Object, testStatistic = testStatistic,
                                      name = name, allMembers = allMembers,
                                      groupMembers = groupMembers)

            .Object@joinFun <- joinFun
            .Object@significant <- sigMembers
            .Object@testStatPar = list(...)

            .Object
          })


setMethod("sigMembers<-", "pC",
          function(object, value) {warning("Assignmanet not allowed"); object})

setMethod("allMembers<-", "pC",
          function(object, value) {warning("Assignmanet not allowed"); object})

setMethod("updateGroup",
          signature(object = "pC", name = "missing", members = "character"),
          function(object, members, parents, sigMembers) {

            ## works for a list with at least one element
            allMembers <- parents[[1]]
            for(s in parents[-1])
              allMembers <- object@joinFun(s, allMembers)

            object@significant <- which(allMembers %in% sigMembers)
            object@allMembers <- allMembers
            object@members <- members

            return(object)
          })

setMethod("updateGroup",
          signature(object = "pC", name = "missing", members = "missing"),
          function(object, parents, sigMembers) {

            allMembers <- parents[[1]]
            for(s in parents[-1])
              allMembers <- object@joinFun(s, allMembers)

            object@significant <- which(allMembers %in% sigMembers)
            object@allMembers <- allMembers

            return(object)
          })




######################################################################
######################################################################

#################### leaCount methods ####################

#################### constructor ####################
setMethod("initialize", "leaCount",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   elim = integer(),
                   depth = 2,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers,
                                      sigMembers, elim = elim)
            .Object@depth <- as.integer(depth)
            .Object@testStatPar = list(...)

            .Object
          })

if(!isGeneric("depth<-"))
  setGeneric("depth<-", function(object, value) standardGeneric("depth<-"))

setMethod("depth<-", "leaCount",
          function(object, value)  {object@depth <- as.integer(value); object})


if(!isGeneric("depth"))
  setGeneric("depth", function(object) standardGeneric("depth"))

setMethod("depth", "leaCount", function(object) object@depth)



#################### elimScore methods ####################

#################### constructor ####################
setMethod("initialize", "leaScore",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   score = numeric(),
                   scoreOrder = "increasing",
                   elim = integer(),
                   depth = 2,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers, score,
                                      scoreOrder, elim = elim)
            .Object@depth <- as.integer(depth)
            .Object@testStatPar = list(...)

            .Object
          })

setMethod("depth<-", "leaScore",
          function(object, value)  {object@depth <- as.integer(value); object})

setMethod("depth", "leaScore", function(object) object@depth)



#################### elimExpr methods ####################

#################### constructor ####################
setMethod("initialize", "leaExpr",
          function(.Object,
                   testStatistic,
                   name = character(),
                   groupMembers = character(),
                   exprDat,
                   pType = factor(),
                   elim = integer(),
                   depth = 2,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      groupMembers, exprDat, pType,
                                      elim = elim)
            .Object@depth <- as.integer(depth)
            .Object@testStatPar = list(...)

            .Object
          })

setMethod("depth<-", "leaExpr",
          function(object, value)  {object@depth <- as.integer(value); object})

setMethod("depth", "leaExpr", function(object) object@depth)
