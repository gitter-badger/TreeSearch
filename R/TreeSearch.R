#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#'  
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
#' @template datasetParam
#' @param outgroup a vector listing the taxa in the outgroup;
#' @param concavity concavity constant for implied weighting (not currently implemented!); 
#' @param method rearrangements to perform; one of \kbd{NNI}, \kbd{SPR}, or \kbd{TBR};
#' @param maxIter the maximum number of iterations to perform before abandoning the search;
#' @param maxHits the maximum times to hit the best pscore before abandoning the search;
#' @param forestSize the maximum number of trees to return - useful in concert with \code{\link{consensus}};
#' @param cluster a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;
#' @param verbosity higher values provide more verbose user feedback in stdout;
#' @param \dots other arguments to pass to subsequent functions.
#' 
#' @return{
#' This function returns a tree, with an attribute \code{pscore} conveying its parsimony score.
#' Note that the parsimony score will be inherited from the tree's attributes, which is only valid if it 
#' was generated using the same \code{data} that is passed here.
#' }
#' @author Martin R. Smith
#'
#' @seealso
#' \itemize{
#' \item \code{\link{Fitch}}, calculates parsimony score;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#' \item \code{\link{SectorialSearch}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{Ratchet}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' data('SigSut')
#' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
#' njtree <- Root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
#'
#' \dontrun{
#' TreeSearch(njtree, SigSut.phy, outgroup, maxIter=20, method='NNI')
#' TreeSearch(njtree, SigSut.phy, outgroup, maxIter=20, method='SPR')
#' TreeSearch(njtree, SigSut.phy, outgroup, maxIter=20, method='TBR')}
#' 
#' @keywords  tree 
#' 
#' @export
TreeSearch <- function 
(tree, dataset, method='NNI', maxIter=100, maxHits=20, forestSize=1, cluster=NULL, 
 verbosity=1, ...) {
  # Initialize morphy object
  if (class(dataset) != 'phyDat') stop ("dataset must be of class phyDat, not ", class(dataset))
  tree <- RenumberTips(tree, names(dataset))
  ret <- DoTreeSearch(tree, dataset, method, maxIter, maxHits, forestSize, cluster, 
                      verbosity, ...)
  return (ret)
}

#' DoTreeSearch
#'
#' Performs a tree search
#' 
#' Does the hard work of searching for a most parsimonious tree.
#' End-users are expected to access this function through its wrapper, TreeSearch
#' It is also called directly by Ratchet and Sectorial functions
#'
#' @template labelledTreeParam
#' @param data dataset in the format expected by TreeScorer
#' @param TreeScorer function to generate optimality score; defaults to \code{\link{FitchScore}}
#' @param Rearrange function to be used for tree rearrangements: probably \link{\code{NNI}},
#'        \link{\code{SPR}} or \link{\code{TBR}}
#' @param maxIter Maximum iterations
#' @param maxHits stop search after finding optimal score \code{maxHits} times
#' @param forestSize number of trees to store in memory
#' @param track Verbosity of reporting
#'
#' @return a tree of class \code{phylo} with attributes "hits" (number of times hit) and "pscore"
#'         (score given by TreeScorer)
#'
#' @author Martin R. Smith
#' 
#' @keywords internal
#' @export
DoTreeSearch <- function (tree, data, TreeScorer = FitchScore, Rearrange = NNI,
                        maxIter = 100, maxHits = 20, forestSize = 1,
                        cluster = NULL, track = 1, ...) {
  tree$edge.length <- NULL # Edge lengths are not supported
  attr(tree, 'hits') <- 1
  if (exists("forestSize") && length(forestSize) && forestSize > 1) {
    forest <- empty.forest <- vector('list', forestSize)
    forest[[1]] <- tree
  } else {
    forestSize <- 1 
  }
  if (is.null(attr(tree, 'score'))) attr(tree, 'score') <- TreeScorer(tree, data)
  bestScore <- attr(tree, 'score')
  if (track > 0) cat("\n  - Performing", method, "search.  Initial score:", bestScore)
  returnSingle <- !(forestSize > 1)
  
  for (iter in 1:maxIter) {
    trees <- RearrangeTree(tree, data, Rearrange, TreeScorer, minScore=bestScore,
                           returnSingle=returnSingle, iter=iter, cluster=cluster, track=track)
    iterScore <- attr(trees, 'score')
    if (length(forestSize) && forestSize > 1) {
      hits <- attr(trees, 'hits')
      if (iterScore == bestScore) {
        forest[(hits-length(trees)+1L):hits] <- trees
        tree <- sample(forest[1:hits], 1)[[1]]
        attr(tree, 'score') <- iterScore
        attr(tree, 'hits') <- hits
      } else if (iterScore < bestScore) {
        bestScore <- iterScore
        forest <- empty.forest
        forest[1:hits] <- trees
        tree <- sample(trees, 1)[[1]]
        attr(tree, 'score') <- iterScore
        attr(tree, 'hits') <- hits
      }      
    } else {
      if (iterScore <= bestScore) {
        bestScore <- iterScore
        tree <- trees
      }
    }
    if (attr(trees, 'hits') >= maxHits) break
  }
  if (track > 0) cat("\n  - Final score", attr(tree, 'score'), "found", attr(tree, 'hits'), "times after", iter, method, "iterations\n")  
  if (forestSize > 1) {
    if (hits < forestSize) forest <- forest[-((hits+1):forestSize)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- bestScore
    return(unique(forest))
  } else {
    return(tree)
  }
}