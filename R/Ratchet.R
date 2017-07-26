#' Parsimony Ratchet
#'
#' \code{Ratchet} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#'
#' @template treeParam 
#' @param data a dataset in the format required by TreeScorer
#' @template concavityParam
#' @param returnAll Set to \code{TRUE} to report all MPTs encountered during the search, perhaps to analyze consensus
#' @param outgroup a vector specifying all tips in the outgroup; if unspecified then identical trees with different roots will be considered unique;
#' @param maxit   maximum ratchet iterations to perform;
#' @param maxIter maximum rearrangements to perform on each bootstrap or ratchet iteration;
#' @param maxHits maximum times to hit best score before terminating a tree search within a pratchet iteration;
#' @param k stop when k ratchet iterations have found the same best score;
#' @param verbosity larger numbers provides more verbose feedback to the user;
#' @param rearrangements method(s) to use when rearranging trees: 
#'        a vector containing any combination of the strings "NNI", "SPR" or "TBR";
#' @param \dots other arguments to pass to subsequent functions.
#' 
#' @return This function returns a tree modified by parsimony ratchet iteration, retaining the position of the root.
#'
#' @references Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.} Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}
#'
#' @author Martin R. Smith
#' 
#' Adapted from \code{\link[phangorn]{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
#' 
#' @seealso \code{\link{Ratchet}}
#' @seealso \code{\link{TreeSearch}}
#' @seealso \code{\link{SectorialSearch}}
#' 
#' @examples Ratchet(RandomTree(Lobo.phy), SigSut.phy, outgroup='Cricocosmia')
#' 
#' @keywords  tree 
#' @export
## TODO use Rooted NNI / SPR / TBR 
Ratchet <- function (tree, data, TreeScorer=FitchScore, returnAll=FALSE, outgroup=NULL, 
                      ratchIter=100, searchIter=5000, searchHits=40, ratchHits=10, track=0, 
                      rearrangements="NNI", suboptimal=1e-08, ...) {
  epsilon <- 1e-08
  if (is.null(attr(tree, "score"))) attr(tree, "score") <- TreeScorer(tree, data)
  bestScore <- attr(tree, "score")
  if (track >= 0) cat("\n* Initial score:", bestScore)
  if (returnAll) {
    null.forest <- vector('list', ratchIter)
    forest <- null.forest
    forest.scores <- rep(NA, ratchIter)
  }

  bestScoreHits <- 0
  iterationsCompleted <- 0
  for (i in 1:ratchIter) {
    if (track >= 0) cat ("\n - Running NNI on bootstrapped dataset. ")
    bootstrapTree <- BootstrapTree(phy=tree, x=data, maxIter=searchIter, maxHits=searchHits,
                        TreeScorer=TreeScorer, track=track-1, ...)
    
    if (track >= 0) cat ("\n - Running", ifelse(is.null(rearrangements), "NNI", rearrangements), "from new candidate tree:")
    if (rearrangements == "TBR") {
      candidate <- DoTreeSearch(bootstrapTree, data, TreeScorer=TreeScorer, method='TBR',
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
      candidate <- DoTreeSearch(candidate, data, TreeScorer=TreeScorer, method='SPR', 
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
      candidate <- DoTreeSearch(candidate, data, TreeScorer=TreeScorer, method='NNI', 
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
    } else if (rearrangements == "TBR only") {  
      candidate <- DoTreeSearch(bootstrapTree, data, TreeScorer=TreeScorer, method='TBR', 
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
    } else if (rearrangements == "SPR") {       
      candidate <- DoTreeSearch(bootstrapTree, data, TreeScorer=TreeScorer, method='SPR', 
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
      candidate <- DoTreeSearch(candidate, data, TreeScorer=TreeScorer, method='NNI', 
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
    } else if (rearrangements == "SPR only") {  
      candidate <- DoTreeSearch(bootstrapTree,    data, TreeScorer=TreeScorer, method='SPR', 
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
    } else {  
      candidate <- DoTreeSearch(bootstrapTree,    data, TreeScorer=TreeScorer, method='NNI', 
                                track=track, maxIter=searchIter, maxHits=searchHits, ...)
    }
    candScore <- attr(candidate, 'score')
    if ((candScore + epsilon) < bestScore) {
      # New 'best' tree
      if (returnAll) {
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup) # TODO Improve this function
        forest.scores[i] <- candScore
      }
      tree <- candidate
      bestScore <- candScore
      bestScoreHits <- 1
    } else if (bestScore + epsilon > candScore) { # i.e. best == cand, allowing for floating point error
      bestScoreHits <- bestScoreHits + 1
      tree <- candidate
      if (returnAll) {
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
        forest.scores[i] <- candScore
      }
    } else if (candScore < (bestScore + suboptimal) && returnAll) {
      forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
      forest.scores[i] <- candScore
    }
    if (track >= 0) cat("\n* Best score after", i, "/", ratchIter, "pratchet iterations:", bestScore, "( hit", bestScoreHits, "/", ratchHits, ")")
    if (bestScoreHits >= ratchHits) {
      iterationsCompleted <- i
      break()
    }
  } # end for
  if (iterationsCompleted == 0) iterationsCompleted <- ratchIter
  if (track >= 0) cat ("\nCompleted parsimony ratchet after", iterationsCompleted, "iterations with score", bestScore, "\n")
   
  if (returnAll) {
    keepers <- !is.na(forest.scores) & forest.scores < bestScore + suboptimal
    forest.scores <- forest.scores[keepers]
    forest <- forest[keepers]
    if (length(forest) > 1) {
      class(forest) <- 'multiPhylo'
      ret <- unique(forest)
    } else if (length(forest) == 1) {
      class(forest) <- 'phylo'
      ret <- forest
    } else {
      stop('No trees!? Is suboptimal set to a sensible (positive) value?')
    }
    scores.unique <- vapply(ret, attr, double(1), 'score')
    cat('Found', sum(scores.unique == min(scores.unique)), 'unique MPTs and', length(ret) - sum(scores.unique == min(scores.unique)), 'suboptimal trees.\n')
    if (is.null(outgroup)) warning('"outgroup" not specified, so some "unique" trees may have same topology but distinct roots.')
  } else {
    ret <- tree
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}