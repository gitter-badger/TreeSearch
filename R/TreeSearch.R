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
    message("  - Performing tree search.  Initial score: ", 
            signif(bestScore, 7))
  }
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) {
    if (verbosity > 0L) {
      message("  - Aborting tree search as tree score ", signif(bestScore, 6), 
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
                signif(scoreThisIteration, 6), " doesn't beat ", 
                signif(bestScore, 6))
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
    message("  - Final score ", signif(bestScore, 7), " found ", hits, " times after ",
            iter, " rearrangements.", if (verbosity > 1L) '\n' else '')
  }
  

  edgeList[3:4] <- c(bestScore, hits)
    
  # Return:
  edgeList
}

#' Shuffle array by third dimension
#' @param array three-dimensional array
#' @author Martin R. Smith
#' @keywords internal
#' @export
ShuffleArray <- function (array) array[, , sample(seq_len(dim(array)[3])), 
                                       drop = FALSE]

#' Bind slices of two arrays
#' @param array1,array2 Three-dimensional arrays
#' @return The two arrays, bound along dimension three
#' @keywords internal
#' @export
BindArrays <- function(array1, array2, dim1 = dim(array1), sliceDim = dim1[1:2],
                       slices1 = dim1[3], slices2 = dim(array2)[3]) {
  array(c(array1, array2), c(sliceDim, slices1 + slices2))
}

#' Write log message
Report <- function (verbosity, level, ..., appendLF = TRUE, appendPrefix = TRUE) {
  if (verbosity > level) { # TODO: Confirm that this really does look for 
    # verbosity in the parent enviroment as it ought to, rather than starting 
    # in the global environment
    prefix <- if (appendPrefix) switch(as.character(level),
                     '0' = '',
                     '1' = ' - ',
                     '2' = ' . ',
                     '3' = '   - ',
                     '4' = '     > ',
                     '     . ') else ''
    message(prefix, ..., appendLF = appendLF)
  }
}

#' @param maxQueue Integer specifying maximum number of candidate trees to 
#' queue for analysis.  Higher numbers use more memory but ensure a more 
#' comprehensive search.
#' @param proposalLimit Integer.  If not `NULL`, `ProposedMoves` will
#' be terminated once `proposalLimit` moves have been proposed.  
#' (`proposalLimit` is sent as the `sampleSize` parameter to `ProposedMoves`).
#' If none of these moves results in a better tree, `proposalLimit` will be
#' increased and more moves requested from `ProposedMoves`.
EdgeMatrixSearch <- function (edgeMatrix, dataset,
                              TreeScorer = MorphyLength,
                              ProposedMoves = RootedTBRSwapAll,
                              maxHits = 40L, maxQueue = 1e06,
                              proposalLimit = 20L,
                              verbosity = 1L, ...) {
  epsilon <- 1e-07
  limitProposals <- !is.null(proposalLimit)
  lastLimit <- 0L
  bestScore <- TreeScorer(edgeMatrix[, 1], edgeMatrix[, 2], dataset, ...)
  
  Report(verbosity, 0L,
         "Performing tree search.  Initial score: ", signif(bestScore, 7))
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
    
    priorHits <- duplicated(BindArrays(
      candidates, hits[, , actualHits, drop = FALSE],
      sliceDim = dimCandidates[1:2], slices1 = nCandidates, slices2 = nHits),
      MARGIN = 3L, fromLast = TRUE)[seq_len(nCandidates)]
    candidates[, , !priorHits, drop = FALSE]
  }
  
  NewCandidates <- function (edgeMatrix) {
    Report(verbosity, 4L, 'Requesting ', proposalLimit, 
           if (!is.null(proposalLimit)) ' ', 'moves... ', appendLF = FALSE)
    
    candidates <- ProposedMoves(edgeMatrix[, 1], edgeMatrix[, 2], nEdge,
                                sampleSize = proposalLimit)
    
    Report(verbosity, 5L, dim(candidates)[3], ' moves proposed, ', 
           appendPrefix = FALSE, appendLF = FALSE)
    
    candidates <- ShuffleArray(NotHitAlready(candidates))
    
    Report(verbosity, 4L, dim(candidates)[3], ' trees returned.',
           appendPrefix = FALSE)
    
    candidates
  }
  
  candidates <- NewCandidates(edgeMatrix)
  lastProposal <- dim(candidates)[3]
  newIteration <- TRUE
  
  while (nHits < maxHits) {
    if (limitProposals && !newIteration) {
      lastLimit <- proposalLimit
      proposalLimit <- lastProposal + lastProposal
      
      Report(verbosity, 4L, 'Increasing proposal limit from ', lastLimit, ' to ',
                proposalLimit)
      
      candidates <- NewCandidates(hits[, , 1])
      if (dim(candidates)[3] < proposalLimit) {
        # We've reached the limit
        Report(verbosity, 4L, 'Proposal limit too high. Removing limit.')
        proposalLimit <- NULL
        limitProposals <- FALSE
      }
      candidates <- candidates[, , -seq_len(lastProposal), drop = FALSE]
      nCandidates <- dim(candidates)[3]
      lastProposal <- nCandidates
    } else {
      nCandidates <- dim(candidates)[3]
    }
    newIteration <- FALSE
    
    if (nCandidates == 0) {
      Report(verbosity, 2L, "No unvisited candidate trees.")
      break
    }
    if (nCandidates > maxQueue) {
      Report(verbosity, 1L, 'Trimming overflowing queue to maxQueue = ', maxQueue)
      candidates <- candidates[, , seq_len(maxQueue)]
    }
    
    Report(verbosity, 4L, nCandidates, ' unvisited trees in queue.')
    
    for (i in seq_len(nCandidates)) {
      candidateScore <- TreeScorer(candidates[, 1, i],candidates[, 2, i],
                                   dataset, ...)
      
      if (candidateScore < bestScore + epsilon) {
        if (candidateScore + epsilon < bestScore) {
          
          bestScore <- candidateScore
          Report(verbosity, 1L, 'New best score ', signif(bestScore), ' on candidate ',
                    i, '/', nCandidates, '.')
          
          if (limitProposals) {
            Report(verbosity, 4L, 'Decreasing proposal limit from ',
                      proposalLimit, appendLF = FALSE)
            
            proposalLimit <- ceiling(proposalLimit / 2L)
            
            Report(verbosity, 4L, ' to ', proposalLimit, appendPrefix = FALSE)
            
            lastProposal <- dim(candidates)[3]
            newIteration <- TRUE
          }
          
          hits <- nothingHit
          hits[, , 1] <- candidates[, , i]
          nHits <- 1L
          candidates <- NewCandidates(hits[, , 1])
          
          i <- 0 # In case we found a hit on the last candidate
          
          break
          
        } else {
          
          nHits <- nHits + 1L
          hits[, , nHits] <- candidates[, , i]
          
          Report(verbosity, 2L, 'Score ', signif(candidateScore, 6), ' hit ',
                    nHits, ' times on candidate ', i, '/', nCandidates, '.')
          
          if (nHits >= maxHits) {
            Report(verbosity, 2L, "Reached maximum hits (maxHits = ", maxHits, ").")
            
            candidates <- array(dim = c(0, 0, 0))
            break
          }
          
          candidates <- BindArrays(NewCandidates(candidates[, , i]),
                                   candidates[, , -seq_len(i), drop = FALSE])
          newIteration <- TRUE
          i <- 0 # In case we found a hit on the last candidate
          break
        }
      } else {
        if (verbosity > 4L) {
          if (verbosity > 8L) { 
            Report(verbosity, 4L, 'Candidate score ', signif(candidateScore, 6), ' > ',
                    signif(bestScore, 6))
          } else {
            if (i %% 20 == 0) message('     ...')
          }
        }
      }
    }
    
    if (i == nCandidates) {
      if (!limitProposals) {
        Report(verbosity, 1L, "All ", nHits, " best trees are locally optimal.")
        break
      } else {
        Report(verbosity, 3L, 'Completed batch of ', nCandidates, ' candidate trees.  ',
               if (newIteration) 'Getting next candidates from queue.' else 
                 'Generating more.')
      }
    }
    
  }
  
  if (verbosity > 0L) {
    message("  - Final score ", signif(bestScore, 7), " found ",
            nHits, " times.", if (verbosity > 1L) '\n' else '')
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


#' @param tree Object of class phylo, corresponding to the starting tree 
#' from which edges were originally obtained.
#' @author Martin R. Smith
#' @keywords internal
#' @export
EdgesToForest <- function (edges, tree) {
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
  ret
}

#' Find local optimum
#' 
#' @inheritParams TreeSearch
#' @author Martin R. Smith
#' @export
FindPeak <- function (tree, dataset, 
                      InitializeData = PhyDat2Morphy,
                      CleanUpData    = UnloadMorphy,
                      TreeScorer     = MorphyLength,
                      ProposedMoves  = RootedTBRSwapAll,
                      maxQueue = 1e06, maxHits = 40L, verbosity = 1L, ...) {
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) {
    stop("tree must be bifurcating; try rooting with ape::root")
  }
  tree <- RenumberTips(tree, names(dataset))
  edgeMatrix <- RenumberTreeStrict(tree$edge[, 1], tree$edge[, 2])
  
  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))
  
  edges <- EdgeMatrixSearch(edgeMatrix, initializedData, TreeScorer=TreeScorer,
                           ProposedMoves = ProposedMoves, maxHits = maxHits,
                           maxQueue = maxQueue,
                           verbosity = verbosity, ...)
  
  # Return:
  EdgesToForest(edges, tree)
}
