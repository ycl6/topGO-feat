#################### topGOdata class ####################
## 
## The node atributes are environments containing the genes/probes
## annotated to the respective node
##
## If genes is a numeric vector than this should represent the gene's
## score. If it is factor it should discriminate the genes in
## interesting genes and the rest

## TODO: it will be a good idea to replace the allGenes and allScore with an
## ExpressionSet class. In this way we can use tests like global test, globalAncova....
##  -- ALL variables sarting with . are just for internal class usage (private)


setClass("topGOdata",
         representation = representation(
           ## some description of the data/experiment
           description = "character",
           ## which of the ontology to be used: BP, CC, MF, (BC, BM, CM, ...) 
           ontology = "character",
           ## the vector containing the genes/probes used in the experiment
           allGenes = "character",
           ## the vector containing the score or the  of the  genes/probes 
           allScores = "ANY",
           ## if each has a score, geneSelectionFun: function(x, ...)
           geneSelectionFun = "function",
           ## which genes are annotated to the specific ontology
           feasible = "logical",
           ## the GO ontology: graph structure and annotations: leaves to root
           graph = "graphNEL",
           ## nodes in the graph with fewer than "nodeSize" genes are removed
           nodeSize = "integer",
           ## expression matrix                                     #!
           expressionMatrix = "ANY",                                #!
           ## phenotype information (groups that shall be compared) #!
           phenotype = "ANY"))                                      #!



######################## topGOresult class ######################

## probably will add more infos with time here
setClass("topGOresult",
         representation = representation(
           ## some description of the data/experiment
           description = "character",
           ## the p-values of the GO terms (named vector)
           score = "numeric",
           ## which test statistic was used 
           testName = "character",
           ## which algorithm was used 
           algorithm = "character",
           ## stats about the genes ...
           geneData = "ANY"))



######################## groupInfo class ########################
## A virtual class containing basic group (GO term) data:
##     gene names, genes scores, etc...

setClass("groupStats",
         representation = representation(
           ## the name of the group (GO ID for example)
           name = "character",
           ## the names of all the mebers (gene ID)
           allMembers = "character",
           ## the names of the mebers in the group (gene ID)
           members = "character",
           ## function containg a statistical test takeing the
           ## parameter the object itself 
           testStatistic = "function",
           testStatPar = "list",
           "VIRTUAL"))


#################### classicCount class #########################
## A class that extends the virtual class groupStats by adding 
## a slot representing the significant members.

setClass("classicCount", contains = "groupStats", 
         representation = representation(
           ## the index of which members are significant  
           significant = "integer"))


##################### classicScore class ########################
## A class that extends the virtual class groupStats by adding 
## a slot representing the score of each gene. (used for KS test)

setClass("classicScore", contains = "groupStats",
         representation = representation(
           ## the score for each member, the most important
           ## member has the highest score
           score = "numeric",
           ## scoreOrder = TRUE (decreasing) the max(score) is considering the best score
           ## scoreOrder = FALSE (increasing) the min(score) is considre the best score
           scoreOrder = "logical"))


##################### classicExpr class #########################
## A class that extends the virtual class groupStats. Contains:
##  - a slot representing the gene expression matrix (for all genes).
##  - a slot representing the phenotipic data (can be a matrix or a vector!).
## Used for GlobalTest.

## the expression matrix will not have row names. The row names will be stored in
## allMembers slot of the parent class

setClass("classicExpr", contains = "groupStats",
         representation = representation(
           ## the expression matrix
           eData = "environment", ## we store the expr matrix in this environment for fast class copying 
           ## scoreOrder = TRUE (decreasing) the max(score) is considering the best score
           pType = "factor"))




##################### weight01Count class #########################
## used for elim2 algorithm
setClass("weight01Count", contains = "classicCount", 
         representation = representation(
           ## the index of which of the group members should be removed
           elim = "integer"))


##################### weight01Score class #########################
setClass("weight01Score", contains = "classicScore",
         representation = representation(
           ## the score for each member, the most important
           ## member has the highest score
           elim = "integer"))


##################### weight01Expr class #########################
setClass("weight01Expr", contains = "classicExpr",
         representation = representation(
           elim = "integer"))



####################### elimCount class #########################
setClass("elimCount", contains = "weight01Count", 
         representation = representation(
           cutOff = "numeric"))


####################### elimScore class #########################
setClass("elimScore", contains = "weight01Score",
         representation = representation(
           cutOff = "numeric")) 


####################### elimExpr class ##########################
setClass("elimExpr", contains = "weight01Expr",
         representation = representation(
           cutOff = "numeric"))



###################### weightCount class ########################
setClass("weightCount", contains = "classicCount", 
         representation = representation(
           weights = "numeric",
           sigRatio = "function",
           ## if more sig nodes ar penalized more (should be part of the sigRatio) 
           penalise = "function",  ## penalise = c("more", "less")
           roundFun = "function"))


#################### parentChild class #########################
##
## the parents information is stored in the class' slots
##
## allMembers - character vector containing all genes annotated in the parent(s).
##              If the group has two or more parents, the genes are
##              concatenated in this vector. The separation between
##              parents is given by the splitIndex.
## splitIndex - the possition in the allMembers where a new parent starts
## joinFun - the way the genes in the parents are joinded: union or intersection

setClass("parentChild", contains = "classicCount", 
         representation = representation(
           splitIndex = "integer",
           joinFun = "character"))

## implements same functionality as parentChild but the join operation
## for the parents is done when the class is build (or updated)
setClass("pC", contains = "classicCount", 
         representation = representation(
           joinFun = "function"))




##################### LEA classes #########################
## used for LEA algorithm
## basically the same class as weight01xxx but we need it for the dispatch mechanism 
setClass("leaCount", contains = "weight01Count",
         representation = representation(
           depth = "integer"))

setClass("leaScore", contains = "weight01Score",
         representation = representation(
           depth = "integer"))

setClass("leaExpr", contains = "weight01Expr",
         representation = representation(
           depth = "integer"))

