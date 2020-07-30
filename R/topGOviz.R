## this file contains code for visualizing the GO DAG and other ploting functions ....

showGroupDensity <- function(object, whichGO, ranks = FALSE, rm.one = TRUE) {

  groupMembers <- genesInTerm(object, whichGO)[[1]]
  
  allS <- geneScore(object, use.names = TRUE)
  if(rm.one)
    allS <- allS[allS < 0.99]
    
  xlab <- "Gene' score"
  if(ranks) {
    allS <- rank(allS, ties.method = "random")
    xlab <- "Gene's rank" 
  }
  
  group <- as.integer(names(allS) %in% groupMembers)
  xx <- data.frame(score = allS,
                   group = factor(group,
                     labels = paste(c("complementary", whichGO), "  (", table(group), ")", sep = "")))
  
  return(densityplot( ~ score | group, data = xx, layout = c(1, 2), xlab = xlab))
}



.ps2eps <- function(filename) {
  system(paste('ps2epsi', paste(filename, ".ps", sep=""), sep = ' '), T, F)
  system(paste('epstool --copy --bbox', paste(filename, ".epsi", sep=""),
               paste(filename, ".eps", sep=""), sep = ' '), T, F)
  system(paste('rm -f', paste(filename, ".ps", sep=""), paste(filename, ".epsi", sep=""), sep = ' '), T, T)
}


if(!isGeneric("printGraph"))
  setGeneric("printGraph",
             function(object, result, firstSigNodes, refResult, ...) standardGeneric("printGraph"))


setMethod("printGraph",
          signature(object = "topGOdata", result = "topGOresult",
                    firstSigNodes = "numeric", refResult = "missing"),
          function(object, result, firstSigNodes = 10, fn.prefix = "",
                   useInfo = "def", pdfSW = FALSE) {

            out.fileName <- paste(fn.prefix, algorithm(result), firstSigNodes, useInfo, sep = '_')              
            ## .DOT.FILE.NAME <<- paste(out.fileName, 'dot', sep = '.')
               
            if(pdfSW)
              pdf(file = paste(out.fileName, 'pdf', sep = '.'), width = 10, height = 10)
            else
              postscript(file = paste(out.fileName, 'ps', sep = '.'))
            
            ## plot the graph to the specified device
            par(mai = rep(0, 4))
            gT <- showSigOfNodes(object, score(result), firstSigNodes = firstSigNodes,
                                 swPlot = FALSE, useInfo = useInfo, plotFunction = GOplot)
            plot(gT$complete.dag)
            dev.off()

            ##if(!pdfSW && .Platform$OS.type == "unix")
            ##  .ps2eps(out.fileName)
            
            cat(out.fileName, ' --- no of nodes: ', numNodes(gT$dag), '\n') 
          })


setMethod("printGraph",
          signature(object = "topGOdata", result = "topGOresult",
                    firstSigNodes = "numeric", refResult = "topGOresult"),
          function(object, result, firstSigNodes = 10, refResult,
                   fn.prefix = "", useInfo = "def", pdfSW = FALSE) {

            out.fileName <- paste(fn.prefix, algorithm(result), algorithm(refResult),
                                  firstSigNodes, useInfo, sep = '_')              
            ## .DOT.FILE.NAME <<- paste(out.fileName, 'dot', sep = '.')
               
            if(pdfSW)
              pdf(file = paste(out.fileName, 'pdf', sep = '.'), width = 10, height = 10)
            else
              postscript(file = paste(out.fileName, 'ps', sep = '.'))
            
            ## plot the graph to the specified device
            par(mai = rep(0, 4))
            wN <- names(sort(score(refResult))[1:firstSigNodes])
            gT <- showSigOfNodes(object, score(result), firstSigNodes = firstSigNodes,
                                 wantedNodes = wN, swPlot = FALSE, useInfo = useInfo,
                                 oldSigNodes = score(refResult), plotFunction = GOplot)
            plot(gT$complete.dag)
            dev.off()
            
            ##if(!pdfSW && .Platform$OS.type == "unix") 
            ##  .ps2eps(out.fileName)
            
            cat(out.fileName, ' --- no of nodes: ', numNodes(gT$dag), '\n') 
          })



## this function will plot the GO DAG or parts of it
## sigNodes:     a named vector of terms p-values, the names are the GO terms
## wantedNodes:  the nodes that we want to find, we will plot this nodes with
##               a different color. The vector contains the names pf the nodes
## oldSigNodes:  used to plot the (new) sigNodes in the same collor range
##               as the old ones
## export.to.dot.file: is a global variable given the name of the output .dot file
GOplot <- function(dag, sigNodes, dag.name = 'GO terms', edgeTypes = T,
                   nodeShape.type = c('box', 'circle', 'ellipse', 'plaintext')[3],
                   genNodes = NULL, wantedNodes = NULL, showEdges = T, useFullNames = F,
                   oldSigNodes = NULL, nodeInfo = NULL) {
    
  if(!missing(sigNodes))
    sigNodeInd = TRUE
  else
    sigNodeInd = FALSE
  
  ## we set the global Graphviz attributes
  graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
  graphAttrs$cluster <- NULL

  #graphAttrs$graph$splines <- FALSE
  
  ## set the node shape
  graphAttrs$node$shape <- nodeShape.type

  ## set the fontsize for the nodes labels
  graphAttrs$node$fontsize <- '14'
  #graphAttrs$node$height <- '1.0'
  #graphAttrs$node$width <- '1.5'

  ## set the local attributes lists
  nodeAttrs <- list()
  edgeAttrs <- list()

  ## try to use adaptive node size
  #nodeAttrs$fixedsize[nodes(dag)] <- rep(FALSE, numNodes(dag))
  
  if(is.null(nodeInfo)) {
    nodeInfo <- character(numNodes(dag))
    names(nodeInfo) <- nodes(dag)
  }
  else
    nodeInfo <- paste('\\\n', nodeInfo, sep = '')
  
  ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
  node.names <- nodes(dag)
  if(!useFullNames)
    nodeAttrs$label <- sapply(node.names,
                              function(x) {
                                return(paste(substr(x, 4, nchar(node.names[1])),
                                             nodeInfo[x], sep = ''))
                              })
  else {
    nodeAttrs$label <- paste(node.names, nodeInfo, sep = '')
    names(nodeAttrs$label) <- node.names
  }

  ## we will change the shape and the color of the nodes that generated the dag
  if(!is.null(wantedNodes)) {
    diffNodes <- setdiff(wantedNodes, genNodes)
    if(length(diffNodes) > 0) {
      nodeAttrs$color[diffNodes] <- rep('lightblue', .ln <- length(diffNodes))
      nodeAttrs$shape[diffNodes] <- rep('circle', .ln)
      nodeAttrs$height[diffNodes] <- rep('0.45', .ln)
      ##nodeAttrs$width[diffNodes] <- rep('0.6', .ln)
      ##nodeAttrs$fixedsize[wantedNodes] <- rep(TRUE, .ln)
    }
  }

  ## we will change the shape and the color of the nodes we want back
  if(!is.null(genNodes)) {
    nodeAttrs$color[genNodes] <- rep('lightblue', .ln <- length(genNodes))
    nodeAttrs$shape[genNodes] <- rep('box', .ln)
    #nodeAttrs$fixedsize[genNodes] <- rep(FALSE, .ln)    
  }
  
  ## we will use different fillcolors for the nodes
  if(sigNodeInd) {
    if(!is.null(oldSigNodes)) {
      old.logSigNodes <- log10(sort(oldSigNodes[nodes(dag)]))
      old.range <- range(old.logSigNodes)
      logSigNodes <- log10(sort(sigNodes))
      logSigNodes[logSigNodes < old.range[1]] <- old.range[1]
      logSigNodes[logSigNodes > old.range[2]] <- old.range[2]

      ## debug:  old.range == range(logSigNodes)
      #if(!identical(all.equal(old.range, range(logSigNodes)), TRUE)){
      #  print(old.range)
      #  print(range(logSigNodes))
      #  stop('some stupid error here :)')
      #}
    }
    else
      old.logSigNodes <- logSigNodes <- log10(sort(sigNodes))
    
    
    sigColor <- round(logSigNodes - range(logSigNodes)[1] + 1)
    old.sigColor <- round(old.logSigNodes - range(old.logSigNodes)[1] + 1)


    mm <- max(sigColor, old.sigColor)
    sigColor <- sigColor + (mm - max(sigColor))

    colorMap <- heat.colors(mm)
    nodeAttrs$fillcolor <- unlist(lapply(sigColor, function(x) return(colorMap[x])))
  }
  
  if(!showEdges)
    graphAttrs$edge$color <- 'white'
  else
    ## if we want to differentiate between 'part-of' and 'is-a' edges
    if(edgeTypes)
      ##    0 for a is_a relation,  1 for a part_of relation
      ## edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 'black', 'red')
      edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 'black', 'black')
  

  ##plot(dag, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)

  return(Rgraphviz::agopen(graph = dag, name = dag.name, attrs = graphAttrs,
                nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs))
}





## this function will plot the GO DAG or parts of it
## nodeCounts:  a matrix 2xNR_of_nodes in which for each node, the first
##              row represent the nr of significant in this node
##              and the second row represent the total nr of genes mapped
##              to this node
##
## wantedNodes: the nodes that we want to find, we will plot this nodes with
##              a different color. The vector contains the names pf the nodes

GOplot.counts <- function(dag, wantedNodes, dag.name = 'GO terms',
                          edgeTypes = T, nodeCounts, showEdges = T) {
  
  if(missing(wantedNodes))
    stop('please give the nodes that you are intrested in')
  
  if(missing(nodeCounts))
    stop('We need the nodeCounts.')

  ## we will plot the sig/all genes for all GO terms
  plotSigChart <- function(curPlot, nodeCounts, wantedNodes) {
    buildDrawing <- function(x, col) {
      force(x)
      y <- x * 100 + 1
      function(node, ur, attrs = list(), radConv = 1) {
        nodeCenter <- getNodeCenter(node)
        pieGlyph(y, xpos = getX(nodeCenter), ypos = getY(nodeCenter),
                 radius = getNodeLW(node), col = col)
        drawTxtLabel(txtLabel(node), getX(nodeCenter), getY(nodeCenter))
        
      }
    }

    drawing <- as.list(1:ncol(nodeCounts))
    .wn <- integer(ncol(nodeCounts))
    names(.wn) <- names(drawing) <- colnames(nodeCounts)
    .wn[wantedNodes] <- 1
    
    drawFun <- lapply(drawing,
                      function(x) {
                        if(.wn[x] == 1)
                          col = c('red', 'lightblue')
                        else
                          col = c('yellow', 'lightgreen')
                        
                        buildDrawing(nodeCounts[, x], col)
                      })
                      
    plot(curPlot, drawNode = drawFun)
    
    ## get the DAG root coordinates
    dagRoot <- getGraphRoot(dag, leafs2root = F)
    parentEnv <- environment()
    rootCenter <- NULL
    lapply(AgNode(curPlot),
           function(x) {
             if(name(x) == dagRoot)
               assign('rootCenter', getNodeCenter(x), envir = parentEnv)
           })
    leftMost <- max(getNodeXY(curPlot)$x)

    legend(leftMost / 1.2, getY(rootCenter),
           legend = c("sig", "all", "sig(wanted)", "all(wanted)"),
           fill = c('yellow', 'lightgreen', 'red', 'lightblue'))

  }
  

  ## we set the global Graphviz attributes
  graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
  graphAttrs$cluster <- NULL

  ## set the fontsize for the nodes labels
  graphAttrs$node$fontsize <- '10'
#  graphAttrs$node$height <- '1.0'
 # graphAttrs$node$width <- '1.5'

  ## set the node shape
  graphAttrs$node$shape <- 'ellipse'

  
  ## set the local attributes lists
  nodeAttrs <- list()
  edgeAttrs <- list()

  ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
  #node.names <- nodes(dag)
  #last.char <- nchar(node.names[1])
  #nodeAttrs$label <- sapply(node.names, substr, 4, last.char)
  nodeAttrs$label <- nodes(dag)
  names(nodeAttrs$label) <- nodes(dag)
  
  if(!showEdges)
    graphAttrs$edge$color <- 'white'
  else
    ## if we want to differentiate between 'part-of' and 'is-a' edges
    if(edgeTypes)
      ##    0 for a is_a relation,  1 for a part_of relation
      #edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 'black', 'red')
      edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 'black', 'black')
  
  dagLayout <- Rgraphviz::agopen(graph = dag, name = dag.name, attrs = graphAttrs,
                      nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs)
  
  plotSigChart(dagLayout, nodeCounts, wantedNodes)
}





## putWN     -- the graph is generated with from the firstSigNodes and the
##              wanted Nodes
## putCL     -- we generate the graph from the nodes given by all previous
##              parameters, plus their children. if putCL = 1 than only the 
##              children are added, if putCL = n we get the nodes form the
##              next n levels.
## type      -- used for ploting pie charts
## swPlot    -- if true the graph is ploted, if not no ploting is done.
## useInfo   -- aditional info to be ploted for a node
showSigOfNodes <- function(GOdata, termsP.value, firstSigNodes = 10, reverse = TRUE,
                            sigForAll = TRUE, wantedNodes = NULL, putWN = TRUE,
                           putCL = 0, type = NULL, showEdges = TRUE, swPlot = TRUE,
                           useFullNames = TRUE, oldSigNodes = NULL,
                           useInfo = c('none', 'pval', 'counts', 'def', 'np', 'all')[1],
                           plotFunction = GOplot, .NO.CHAR = 20) {

  requireNamespace('Rgraphviz', quietly = TRUE) || stop('package Rgraphviz is required')

  if(!is.null(firstSigNodes)) 
    sigTerms <- sort(termsP.value)[1:firstSigNodes]
  else
    sigTerms <- numeric(0)
  
  if(putWN && !is.null(wantedNodes))
    baseNodes <- union(names(sigTerms), wantedNodes)
  else
    baseNodes <- names(sigTerms)

  if(length(baseNodes) == 0)
    stop('No nodes were selected')
  
  ## we want to get aditional nodes
  if(putCL) {
    goDAG.r2l <- reverseArch(graph(GOdata))

    for(i in 1:putCL) {
      newNodes <- unique(unlist(adj(goDAG.r2l, baseNodes)))
      baseNodes <- union(newNodes, baseNodes)
    }
  }

  dag <- inducedGraph(graph(GOdata), baseNodes)

  if(reverse)
    dag <- reverseArch(dag)

  termCounts <- termStat(GOdata, nodes(dag))
  
  ## we plot for each node of GO graph the pie plot showing the
  ## difference bettween all genes mapped to it and sig genes mapped to it
  if(!is.null(type)) {
    if(swPlot)
      GOplot.counts(dag, wantedNodes = wantedNodes, nodeCounts = termCounts,
                    showEdges = showEdges)
    return(dag)
  }

  
  pval.info <- function(whichNodes) {
    ret.val <- format.pval(termsP.value[whichNodes], digits = 3, eps = 1e-20)
    names(ret.val) <- whichNodes
    return(ret.val)
  }

  .pval = pval.info(nodes(dag))
  .def = .getTermsDefinition(whichTerms = nodes(dag), ontology(GOdata), numChar = .NO.CHAR)
  .counts = apply(termCounts[, c("Significant", "Annotated")], 1, paste, collapse = " / ")
  ## more infos will be added
  nodeInfo <- switch(useInfo,
                     none = NULL,
                     pval = .pval,
                     def = .def,
                     counts = .counts,
                     np = paste(.def, .pval, sep = '\\\n'),
                     all = paste(.def, .pval, .counts, sep = '\\\n')
                     )
    
  ## we can plot the significance level of all nodes in the dag or for the sigNodes
  if(sigForAll)
    sigNodes <- termsP.value[nodes(dag)]
  else
    sigNodes <- sigTerms

  if(is.null(wantedNodes))
    wantedNodes <- names(sigTerms)

  
  complete.dag <- plotFunction(dag, sigNodes = sigNodes, genNodes = names(sigTerms),
                               wantedNodes = wantedNodes, showEdges = showEdges,
                               useFullNames = useFullNames, oldSigNodes = oldSigNodes,
                               nodeInfo = nodeInfo)
  
  if(swPlot && !is.null(complete.dag))
    plot(complete.dag)
  
  ## we return the obtained dag
  return(list(dag = dag, complete.dag = complete.dag))
}






## this function return a named vector:
## the elements are the edge weights: w(a, b) for a graph
## the names are the edge names in the form: a~b

.getEdgeWeights <- function (graph) {
  
  weightsList <- edgeWeights(graph)
  to <- lapply(weightsList, names)
  from <- nodes(graph)

  if (any(is.na(unlist(to))) || any(is.na(from))) 
    stop("Edge names do not match node names.")

  edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
  edge.weights <- unlist(weightsList)
  names(edge.weights) <- edge.names

  return(edge.weights)
}




## this function is compiling a .dot file from the dag
printDOT <- function(dag, sigNodes = NULL, genNodes = NULL, wantedNodes = NULL,
                     showEdges = T, useFullNames = F, oldSigNodes = NULL,
                     nodeInfo = NULL, export.to.dot.file = "MyGraph.dot") {
  
  ## we set the global Graphviz attributes
  graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
  graphAttrs$cluster <- NULL

  #graphAttrs$graph$splines <- FALSE
  graphAttrs$graph$size <- "6.99,3.99"
  
  ## set the node shape
  graphAttrs$node$shape <- 'ellipse'

  ## set the fontsize for the nodes labels
  graphAttrs$node$fontsize <- '9'
  graphAttrs$edge$fontsize <- '9'
  graphAttrs$node$style <- 'filled'
  #graphAttrs$node$height <- '1.0'
  #graphAttrs$node$width <- '1.5'

  ## set the local attributes lists
  nodeAttrs <- list()
  edgeAttrs <- list()

  ## try to use adaptive node size
  #nodeAttrs$fixedsize[nodes(dag)] <- rep(FALSE, numNodes(dag))
  
  if(is.null(nodeInfo)) {
    nodeInfo <- character(numNodes(dag))
    names(nodeInfo) <- nodes(dag)
  }
  else
    nodeInfo <- paste('\\n<', nodeInfo, '>', sep = '')
  
  ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
  node.names <- nodes(dag)
  if(!useFullNames)
    nodeAttrs$label <- sapply(node.names,
                              function(x) {
                                return(paste(substr(x, 4, nchar(node.names[1])),
                                             nodeInfo[x], sep = ''))
                              })
  else {
    nodeAttrs$label <- paste(node.names, nodeInfo, sep = '')
    names(nodeAttrs$label) <- node.names
  }

  ## we will change the shape and the color of the nodes that generated the dag
  if(!is.null(wantedNodes)) {
    diffNodes <- setdiff(wantedNodes, genNodes)
    if(length(diffNodes) > 0) {
      nodeAttrs$color[diffNodes] <- rep('lightblue', .ln <- length(diffNodes))
      nodeAttrs$shape[diffNodes] <- rep('circle', .ln)
      nodeAttrs$height[diffNodes] <- rep('0.45', .ln)
      ##nodeAttrs$width[diffNodes] <- rep('0.6', .ln)
      ##nodeAttrs$fixedsize[wantedNodes] <- rep(TRUE, .ln)
    }
  }

  ## we will change the shape and the color of the nodes we want back
  if(!is.null(genNodes)) {
    nodeAttrs$color[genNodes] <- rep('lightblue', .ln <- length(genNodes))
    nodeAttrs$shape[genNodes] <- rep('box', .ln)
    #nodeAttrs$fixedsize[genNodes] <- rep(FALSE, .ln)    
  }
  
  ## we will use different fillcolors for the nodes
  if(!is.null(sigNodes)) {
    if(!is.null(oldSigNodes)) {
      old.logSigNodes <- log10(sort(oldSigNodes[nodes(dag)]))
      old.range <- range(old.logSigNodes)
      logSigNodes <- log10(sort(sigNodes))
      logSigNodes[logSigNodes < old.range[1]] <- old.range[1]
      logSigNodes[logSigNodes > old.range[2]] <- old.range[2]
    }
    else
      old.logSigNodes <- logSigNodes <- log10(sort(sigNodes))
    
    sigColor <- round(logSigNodes - range(logSigNodes)[1] + 1)
    old.sigColor <- round(old.logSigNodes - range(old.logSigNodes)[1] + 1)
    
    mm <- max(sigColor, old.sigColor)
    sigColor <- sigColor + (mm - max(sigColor))
    
    colorMap <- heat.colors(mm)
    nodeAttrs$fillcolor <- unlist(lapply(sigColor, function(x) return(colorMap[x])))
  }
  
  if(!showEdges)
    graphAttrs$edge$color <- 'white'
  else
    ## if we want to differentiate between 'part-of' and 'is-a' edges
    ##    0 for a is_a relation,  1 for a part_of relation
    ##  edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 'black', 'red')
    edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 'black', 'black')
  

  toDot(graph = dag, filename = export.to.dot.file,
        attrs = graphAttrs, nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs)

  return(NULL)
}


