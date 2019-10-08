#' @describeIn TreeSearch Tree Search from Edge lists
#' @template edgeListParam
#' @template dataForFunction
#' @author Martin R. Smith
#' @keywords internal
#' @export
EdgeListSearch <- function (edgeList, dataset,
                            TreeScorer = MorphyLength,
                            EdgeSwapper = RootedTBRSwap,
                            maxIter = 100L, maxHits = 20L, 
                            bestScore = NULL, stopAtScore = NULL, 
                            stopAtPeak = FALSE, stopAtPlateau = 0L,
                            verbosity = 1L, ...) {
  epsilon <- 1e-07
  
  if (is.null(bestScore)) {
    if (length(edgeList) < 3L) {
      bestScore <- TreeScorer(edgeList[[1]], edgeList[[2]], dataset, ...)
    } else {
      bestScore <- edgeList[[3]]
    }
  }
  if (verbosity > 0L) {
    message("  - Performing tree search.  Initial score: ", bestScore)
  }
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) {
    if (verbosity > 0L) {
      message("  - Aborting tree search as tree score ", bestScore, 
              " already below target of ", stopAtScore)
    }
    edgeList[[3]] <- bestScore
    return(edgeList)
  }
  
  hits <- 0L
  unimprovedSince <- 0L
  
  for (iter in 1:maxIter) {
    candidateLists <- RearrangeEdges(edgeList[[1]], edgeList[[2]], dataset=dataset, 
                             TreeScorer=TreeScorer, EdgeSwapper=EdgeSwapper, 
                             hits=hits, iter=iter, verbosity=verbosity, ...)
    scoreThisIteration <- candidateLists[[3]]
    hits <- candidateLists[[4]]
  
    if (scoreThisIteration < bestScore + epsilon) {
      if (scoreThisIteration + epsilon < bestScore) unimprovedSince <- -1L
      bestScore <- scoreThisIteration
      edgeList  <- candidateLists
      if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(edgeList)
    } else if (stopAtPeak && scoreThisIteration > bestScore + epsilon) {
      if (verbosity > 1L) {
        message("    ! Iteration ", iter, 
                " - No TBR rearrangement improves score. ",
                scoreThisIteration, " doesn't beat ", bestScore)
      }
      break
    }
    unimprovedSince <- unimprovedSince + 1L
    if (stopAtPlateau > 0L) {
      if (verbosity > 2L && unimprovedSince > 0L) {
        message(" Last improvement ", unimprovedSince, " iterations ago.")
      }
      if (unimprovedSince >= stopAtPlateau) {
        if (verbosity > 1L) message("  - Terminating search, as score has ",
                                    "not improved over past ",
                                    unimprovedSince, " searches.")
        break
      }
    }
  
    if (hits >= maxHits) {
      if (verbosity > 1L) message("  - Terminating search; hit best score ",
                                  hits, " times.")
      break
    }
  }
  if (verbosity > 0L) {
    message("  - Final score ", bestScore, " found ", hits, " times after ",
            iter, " rearrangements.", if (verbosity > 1L) '\n' else '')
  }
  

  edgeList[3:4] <- c(bestScore, hits)
    
  # Return:
  edgeList
}

#' @param followPlateau Logical: if `TRUE`, search all trees one TBR move away
#' from any best-scoring tree; if `FALSE`, stop once a best tree has been found 
#' and no immediate neighbours have a better score.
EdgeMatrixSearch <- function (edgeMatrix, dataset,
                              TreeScorer = MorphyLength,
                              ProposedMoves = RootedTBRSwapAll,
                              maxHits = 40L, followPlateau = TRUE,
                              bestScore=NULL, verbosity=1L, ...) {
  epsilon <- 1e-07
  if (is.null(bestScore)) {
    attrScore <- attr(edgeMatrix, 'score')
    if (is.null(attrScore)) {
      bestScore <- TreeScorer(edgeMatrix[, 1], edgeMatrix[, 2], dataset, ...)
    } else {
      bestScore <- attrScore
    }
  }
  if (verbosity > 0L) {
    message("  - Performing tree search.  Initial score: ", bestScore)
  }
  nHits <- 1L
  nEdge <- dim(edgeMatrix)[1]
  nothingHit <- array(NA, dim = c(nEdge, 2L, maxHits))
  hits <- nothingHit
  hits[, , 1] <- edgeMatrix
  hitToExamine <- 1L
  
  NotHitAlready <- function (candidates) {
    dimCandidates <- dim(candidates)
    nCandidates <- dimCandidates[3]
    actualHits <- !is.na(hits[1, 1, ])
    
    priorHits <- duplicated(array(
      c(candidates, hits[, , actualHits, drop = FALSE]),
      dim = c(dimCandidates[1:2], nCandidates + nHits)),
      MARGIN = 3L, fromLast = TRUE)[seq_len(nCandidates)]
    candidates[, , !priorHits]
  }
  
  while (nHits < maxHits) {
    startTree <- hits[, , hitToExamine]
    candidates <- ProposedMoves(startTree[, 1], startTree[, 2], nEdge)
    candidates <- NotHitAlready(candidates)
    
    nCandidates <- dim(candidates)[3]
    if (verbosity > 3L) {
      message('  - ', nCandidates, ' unvisited trees one TBR move away.')
    }
    stuck <- TRUE
    
    for (i in sample(seq_len(nCandidates))) {
      candidateScore <- TreeScorer(candidates[, 1, i],candidates[, 2, i],
                                   dataset, ...)
      
      if (candidateScore < bestScore + epsilon) {
        stuck <- FALSE
        if (candidateScore + epsilon < bestScore) {
          
          bestScore <- candidateScore
          if (verbosity > 1L) {
            message('  - New best score: ', bestScore)
          }
          
          hits <- nothingHit
          hits[, , 1] <- candidates[, , i]
          nHits <- 1L
          hitToExamine <- 1L
          
          break
          
        } else {
          
          nHits <- nHits + 1L
          hits[, , nHits] <- candidates[, , i]
          
          if (verbosity > 2L) {
            message('   - Score ', candidateScore, ' hit ', nHits, ' times.')
          }
          
          if (nHits >= maxHits) {
            if (verbosity > 2L) {
              message("   - Reached maximum hits (maxHits = ", maxHits, ").")
            }
            break
          }
          
        }
      } else {
        if (verbosity > 4L) {
          message('   - Candidate score ', candidateScore, ' > ',
                  bestScore)
        }
      }
    }
    if (stuck) {
      if (followPlateau) {
        if (hitToExamine < nHits) {
          if (verbosity > 2L) {
            message("  - Result ", hitToExamine, " is locally optimal. ",
                    "Continuing from result ", hitToExamine + 1L, "/",
                    nHits)
          }
          hitToExamine <- hitToExamine + 1L
        } else {
          if (verbosity > 1L) {
            message("  - All ", nHits, " best trees are locally optimal.")
          }
          
          break  
        }
      } else {
        
        if (verbosity > 1L) {
          message("  - Tree is locally optimal. ",
                  "Stopping, as followPlateau = FALSE")
        }
        
        break
      }
    }
  }
  
  if (verbosity > 0L) {
    message("  - Final score ", bestScore, " found ", nHits, " times.",
            if (verbosity > 1L) '\n' else '')
  }
  
  ret <- hits[, , !is.na(hits[1, 1, ]), drop = FALSE]
  attr(ret, 'score') <- bestScore
  
  # Return:
  ret
}

#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#'  
#' For detailed documentation of the TreeSearch package, including full 
#' instructions for loading phylogenetic data into R and initiating and 
#' configuring tree search, see the 
#' [package documentation](https://ms609.github.io/TreeSearch).
#'  
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
#' @template datasetParam
#' @template EdgeSwapperParam
#' @param maxIter the maximum number of iterations to perform before abandoning the search.
#' @param maxHits the maximum times to hit the best pscore before abandoning the search.
#' @template stopAtPeakParam
#' @template stopAtPlateauParam
#'
#' @template InitializeDataParam
#' @template CleanUpDataParam
#' @template treeScorerParam
#'
#' @template verbosityParam
#' @template treeScorerDots
#' 
#' @return {
#' This function returns a tree, with an attribute \code{pscore} conveying its parsimony score.
#' Note that the parsimony score will be inherited from the tree's attributes, which is only valid if it 
#' was generated using the same \code{data} that is passed here.
#' }
#' @author Martin R. Smith
#'
#' @references
#' - \insertRef{SmithTern}{TreeSearch}
#'
#' @seealso
#' \itemize{
#' \item \code{\link{Fitch}}, calculates parsimony score;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#### \item \code{\link{Sectorial}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{Ratchet}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' data('Lobo')
#' njtree <- NJTree(Lobo.phy)
#'
#' \dontrun{
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=NNISwap)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=RootedSPRSwap)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=TBRSwap)}
#' 
#' @keywords  tree 
#' @importFrom Rdpack reprompt
#' @export
TreeSearch <- function (tree, dataset,
                        InitializeData = PhyDat2Morphy,
                        CleanUpData    = UnloadMorphy,
                        TreeScorer     = MorphyLength,
                        EdgeSwapper    = RootedTBRSwap,
                        maxIter = 100L, maxHits = 20L,
                        stopAtPeak = FALSE, stopAtPlateau = 0L,
                        verbosity = 1L, ...) {
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) {
    stop("tree must be bifurcating; try rooting with ape::root")
  }
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- MatrixToList(tree$edge)
  edgeList <- RenumberEdges(edgeList[[1]], edgeList[[2]])

  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))

  edgeList <- EdgeListSearch(edgeList, initializedData, TreeScorer=TreeScorer, 
                             EdgeSwapper = EdgeSwapper, maxIter = maxIter, 
                             maxHits = maxHits, stopAtPeak = stopAtPeak, 
                             stopAtPlateau = stopAtPlateau,
                             verbosity = verbosity, ...)
  
  tree$edge <- ListToMatrix(edgeList)
  attr(tree, 'score') <- edgeList[[3]]
  attr(tree, 'hits') <- edgeList[[4]]
  # Return:
  tree 
}

FindPeak <- function (tree, dataset, 
                         InitializeData = PhyDat2Morphy,
                         CleanUpData    = UnloadMorphy,
                         TreeScorer     = MorphyLength,
                         ProposedMoves  = RootedTBRSwapAll,
                         followPlateau = TRUE,
                         maxHits = 40L, verbosity = 1L, ...) {
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) {
    stop("tree must be bifurcating; try rooting with ape::root")
  }
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- MatrixToList(tree$edge)
  edgeList <- RenumberEdges(edgeList[[1]], edgeList[[2]])
  edgeMatrix <- array(unlist(edgeList), dim = dim(tree$edge))
  
  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))
  
  bestScore <- attr(tree, 'score')
  edges <- EdgeMatrixSearch(edgeMatrix, initializedData, TreeScorer=TreeScorer,
                           ProposedMoves = ProposedMoves, maxHits = maxHits,
                           followPlateau = followPlateau,
                           verbosity = verbosity, ...)
  
  if (dim(edges)[3] > 1L) {
     ret <- structure(apply(edges, 3, function (edge) {
        ret <- tree
        ret$edge <- edge
        ret
      }), class='multiPhylo', score = attr(edges, 'score'))
  } else {
    tree$edge <- edges[, , 1]
    attr(tree, 'score') <- attr(edges, 'score')
    ret <- tree
  }
  
  # Return:
  ret
}