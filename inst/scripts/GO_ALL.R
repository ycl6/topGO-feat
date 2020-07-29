## sample session to exemplify the use of topGO


library(topGO)


## load the data set 
library(ALL)
data(ALL)

## get the library
affyLib <- "hgu95av2"
library(package = affyLib, character.only = TRUE)

## use all the probes that are present on the chip
curr.probes <- intersect(ls(get(paste(affyLib, "ACCNUM", sep = ""))), featureNames(ALL))
ALL <- ALL[curr.probes, ]

edata <- exprs(ALL)
pd <- pData(ALL)




######################################################################
## generate a random list of interesting  genes

## get the list of all genes
geneNames <- featureNames(ALL)
## select (or define) the list of interesting genes
myInterestedGenes <- sample(geneNames, 100)

#str(geneNames)
#str(myInterestedGenes)

## make a indicator vector showing which genes are interesting
geneList <- geneNames %in% myInterestedGenes
sum(geneList)
geneList <- factor(as.integer(geneList))
names(geneList) <- geneNames
str(geneList)

## build the topGOdata class 
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              description = "GO analysis of ALL data WITH MY SELECTED GENES",
              annot = annFUN.hgu,
              affyLib = affyLib)

## display the GOdata object
GOdata




######################################################################
## generate a list of genes

## discriminate B-cell from T-cell
discrGenes <- function(eset = ALL) {
  Tcell.type <- sapply(eset$BT, function(x) return(substr(x, 1, 1) == 'T'))
  return(as.integer(Tcell.type))
}


## the phenotype
y <- discrGenes(pd)

## Differentially expressed genes
geneList <- getPvalues(edata, classlabel = y, alternative = "two.sided")

## function to select the "significant" genes
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
}



######################################################################
## start the GO analysis
## build the topGOdata class 

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSel = topDiffGenes,
              description = "GO analysis of ALL data. Diff. Expre. between B-cell and T-cell",
              annot = annFUN.hgu,
              affyLib = affyLib)

## display the GOdata object
GOdata


#################### Some examples ####################

## description of the experiment
description(GOdata)

## obtain the genes
a <- genes(GOdata)
str(a)
numGenes(GOdata)


## obtain the score of the genes
selGenes <- rownames(edata)[sample(1:nrow(edata), 10)]
gs <- geneScore(GOdata, whichGenes = selGenes)
print(gs)

## if we want a named vector -- this verssion is prefered when whichGenes is used
gs <- geneScore(GOdata, whichGenes = selGenes, use.names = TRUE)
print(gs)

## the list of significant genes
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)

## if we want to update the genes - lets say only the feasible one
.geneList <- geneScore(GOdata, use.names = TRUE)
GOdata ## more available genes
GOdata <- updateGenes(GOdata, .geneList, topDiffGenes)
GOdata ## the available genes are now the feasible genes

## to list the available GO terms (all the nodes in the graph)
str(usedGO(GOdata))

## to list the genes annotated to a set of specified GO terms
sel.terms <- sample(usedGO(GOdata), 10)
ann.genes <- genesInTerm(GOdata, sel.terms)
str(ann.genes)

## if you also need the score for these genes
ann.score <- scoresInTerm(GOdata, sel.terms)
## ann.score <- scoresInTerm(GOdata, sel.terms, use.names = TRUE)
str(ann.score)

## to see the number of annotated genes
num.ann.genes <- countGenesInTerm(GOdata)
str(num.ann.genes)

## to summarise the statistics
termStat(GOdata, sel.terms)




################################################################################



## define the type of test you want to use
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
## run a algorithm that fits the statistic
resultFis <- getSigGroups(GOdata, test.stat)


## ELIM, standar cutOff parameter
#test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
## run a algorithm that fits the statistic
#resultElim.c1<- getSigGroups(GOdata, test.stat)

test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 10)
## run a algorithm that fits the statistic
resultElim <- getSigGroups(GOdata, test.stat)

## for sigRatio use only: ratio, log or 01
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
## run a algorithm that fits the statistic
resultWeight <- getSigGroups(GOdata, test.stat)


## define the type of test you want to use
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
## run a algorithm that fits the statistic
resultKS <- getSigGroups(GOdata, test.stat)

l <- list(classic = score(resultFis),
          elim = score(resultElim),
          weight = score(resultWeight),
          KS = score(resultKS))

##save.image(file = "GO_ALL.RData")

a <- genTable(GOdata, l, orderBy = "weight", ranksOf = "classic", top = 50, numChar = 40, use.levels = FALSE)
## if you want a latex table
if(require(xtable))
  print(xtable(apply(a, 2, as.character)), floating = FALSE)

## to plot a subgraph
showSigOfNodes(GOdata, score(resultFis), firstTerms = 15, useInfo = 'all')


## to print it
##printGraph(GOdata, resultWeight, 15)

printGraph(GOdata, resultWeight, firstSigNodes = 10, fn.prefix = "ALL_BT")

## print the diff between to methods
printGraph(GOdata, resultWeight, firstSigNodes = 15, resultFis, fn.prefix = "ALL_BT", useInfo = "all")
printGraph(GOdata, resultElim, firstSigNodes = 15, resultFis, fn.prefix = "ALL_BT", useInfo = "all")

