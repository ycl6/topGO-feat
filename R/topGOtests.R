## put here all tests statistic

if(!isGeneric("GOFisherTest"))
  setGeneric("GOFisherTest", function(object) standardGeneric("GOFisherTest"))

setMethod("GOFisherTest", "classicCount",
          function(object) {

            contMat <- contTable(object)

            if(all(contMat == 0))
              p.value <- 1
            else
              p.value <- fisher.test(contMat, alternative = 'greater')$p.value

            return(p.value)
          })



if(!isGeneric("GOKSTest"))
  setGeneric("GOKSTest", function(object) standardGeneric("GOKSTest"))

setMethod("GOKSTest", "classicScore",
          function(object) {

            N <- numAllMembers(object)
            na <- numMembers(object)

            ## if the group is empty ... (should not happen, but you never know!)
            if(na == 0 || na == N)
              return(1)

            x.a <- rankMembers(object)
            ## x.b <- setdiff(1:N, x.a)

            return(ks.test(x.a, seq_len(N)[-x.a], alternative = "greater")$p.value)
          })



if(!isGeneric("GOtTest"))
  setGeneric("GOtTest", function(object) standardGeneric("GOtTest"))

setMethod("GOtTest", "classicScore",
          function(object) {

            ## we should try to do this smarter!!!
            var.equal <- FALSE
            ## check if there are any parameters set
            if(length(testStatPar(object)) != 0)
              if("var.equal" %in% names(testStatPar(object)))
                var.equal <- testStatPar(object)[["var.equal"]]

            nG <- numMembers(object)
            nNotG <- numAllMembers(object) - nG

            ## if the group/complementary is empty
            if(nNotG == 0 || nG == 0)
              return(1)

            x.G <- membersScore(object)
            x.NotG <- allScore(object, TRUE)[setdiff(allMembers(object), members(object))]

            aa <- scoreOrder(object)
            ## alternative - TRUE if the max(score) is the best score), FALSE if the min(score)...
            if(nNotG == 1) {
              .stat <- ifelse(aa, sum(x.G >= x.NotG), sum(x.G <= x.NotG))
              p.value <- (.stat + 1) / (nG + nNotG + 1)
            }
            else if(nG == 1) {
              .stat <- ifelse(aa, sum(x.NotG >= x.G), sum(x.NotG <= x.G))
              p.value <- (.stat + 1) / (nG + nNotG + 1)
            }
            else
              p.value <- t.test(x = x.G, y = x.NotG, var.equal = var.equal,
                                alternative = ifelse(aa, "greater", "less"))$p.value

            return(p.value)
          })



if(!isGeneric("GOglobalTest"))
  setGeneric("GOglobalTest", function(object) standardGeneric("GOglobalTest"))

setMethod("GOglobalTest", "classicExpr",
          function(object) {

            if(numMembers(object) == 0)
              return(1)

            requireNamespace('globaltest', quietly = TRUE)
	    return(globaltest::p.value(globaltest::gt(X = membersExpr(object), Y = pType(object))))
          })


permSumStats <- function(object, N) {

  ## the gene scores
  scoreVec <- geneScore(object)
  sizeLookUp <- integer(length(scoreVec))

  ## the available GOs and their size
  goSize <- sort(unique(countGenesInTerm(object)))

  ## remove the GO sizes which are 60% or more of the max size
  goSize <- goSize[goSize / max(goSize) < .6]
  maxSample <- goSize[length(goSize)]

  ## update the look-up table with the feasible sizes
  sizeLookUp[goSize] <- 1:length(goSize)

  assign(".PERMSUM.LOOKUP", sizeLookUp, envir = as.environment(pos = 1L))
  assign(".PERMSUM.MAT", sapply(1:N, function(x) cumsum(sample(scoreVec, maxSample))[goSize]), envir = as.environment(pos = 1L))
}

## same as above just doesn't remove any term
permSumStats.all <- function(object, N) {

  ## the gene scores
  scoreVec <- geneScore(object)
  sizeLookUp <- integer(length(scoreVec))

  ## the available GOs and their size
  goSize <- sort(unique(countGenesInTerm(object)))

  ## update the look-up table with the feasible sizes
  sizeLookUp[goSize] <- 1:length(goSize)

  assign(".PERMSUM.LOOKUP", sizeLookUp, envir = as.environment(pos = 1L))
  assign(".PERMSUM.MAT", sapply(1:N, function(x) cumsum(sample(scoreVec))[goSize]), envir = as.environment(pos = 1L))
}


## Category test based on score sums
if(!isGeneric("GOSumTest"))
  setGeneric("GOSumTest", function(object) standardGeneric("GOSumTest"))

setMethod("GOSumTest", "classicScore",
          function(object) {

            nG <- numMembers(object)
            if(nG == 0) return(1)
            N <- ncol(.PERMSUM.MAT)

            ## compute the group sum (statistic)
            sObserved <- sum(membersScore(object))

            ## check to see if we have a precomputed the permutations for this group size
            isPermuted <- .PERMSUM.LOOKUP[nG]
            if(isPermuted > 0)
              ## obtain the vector of already permuted values and compute the p-value
              return((sum(.PERMSUM.MAT[isPermuted, ] >= sObserved) + 1) / (N + 1))

            ## at this point we don't have permuations for the current group size
            numScore <- numAllMembers(object)  ##  length(scoreVec)
            scoreVec <- allScore(object, FALSE)

            ##cat("\rComputing", N, "permutations for group size:", nG, "...\n")
            ##flush.console()

            if(nG > numScore / 2) {
              ## we simply compute the sum statistic
              permSums <- sapply(1:N, function(x) sum(sample(scoreVec, nG)))
            }
            else { ## we try to do a bit of optimisation using the cumsum function
              chunkIndex <- seq.int(from = nG, to = numScore, by = nG)
              cycles <- N %/% (length(chunkIndex) - 1) + 1

              permSums <-  unlist(lapply(1:cycles, function(y) {
                x <- cumsum(scoreVec[.Internal(sample(numScore, numScore, FALSE, NULL))])[chunkIndex]
                return(x[-1] - x[-length(x)])
              }), use.names = FALSE)

              length(permSums) <- N
            }

            return((sum(permSums >= sObserved) + 1) / (N + 1))
          })

if(!isGeneric("GOKSTiesTest"))
  setGeneric("GOKSTiesTest", function(object) standardGeneric("GOKSTiesTest"))

setMethod("GOKSTiesTest", "classicScore",
          function(object) {

            N <- numAllMembers(object)
            na <- numMembers(object)

            ## if the group is empty ... (should not happen, but you never know!)
            if(na == 0 || na == N)
              return(1)

            ## get the test parameters ... if existent
            numPerm <- testStatPar(object)[["numPerm"]]
            FUN <- testStatPar(object)[["FUN"]]
            if(is.null(numPerm))
              numPerm <- 50L
            if(is.null(FUN))
              FUN <- "max"

            ## get all the scores and take care of their order
            a.s <- allScore(object)
            ## get the members index
            mem.index <- rankMembers(object)
            seqN <- seq_len(N)
            resP <- numeric(numPerm)
            ## for each permutation we get random ranks for the ties
            for(i in 1:numPerm) {
              x.a <- sort.list(order(a.s, sample(N), decreasing = scoreOrder(object)))[mem.index]
              resP[i] <- ks.test(x.a, seqN[-x.a], alternative = "greater")$p.value
            }

            ## aggregate the resulting p-values and return
            return(get(FUN)(resP))
          })

