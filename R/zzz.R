##.onLoad <- function(lib, pkg) {
##  require(methods)
##}

.onAttach <- function(lib, pkg) {
  ## some preprocessing
  where <- match(paste("package:", pkg, sep=""), search())
  groupGOTerms(where)
}
