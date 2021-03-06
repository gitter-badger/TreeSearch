#' Nearest Neighbour Interchange (NNI)
#' 
#' Performs a single iteration of the nearest-neighbour interchange algorithm.
#' Based on the corresponding \code{phangorn} function, but re-coded to improve speed.
#' 
#' Branch lengths are not supported.
#' 
#' All nodes in a tree must be bifurcating; [ape::collapse.singles] and
#' [ape::multi2di] may help.
#' 
#' @template treeParam
#' @template edgeToBreakParam
#' 
#' @return Returns a tree with class \code{phylo} (if \code{returnAll = FALSE}) or 
#'         a set of trees, with class \code{multiPhylo} (if \code{returnAll = TRUE}).
#'
#' @references
#' The algorithm is summarized in
#'  \insertRef{Felsenstein2004}{TreeSearch}
#' 
#' @author Martin R. Smith
#' @family tree rearrangement functions
#' 
#' @examples
#' tree <- ape::rtree(20, br = NULL)
#' NNI(tree)
#' NNI(tree, edgeToBreak = -1)
#'
#' @export
NNI <- function (tree, edgeToBreak=NULL) {
  edge <- tree$edge
  parent <- edge[, 1]
  StopUnlessBifurcating(parent)
  if (!is.null(edgeToBreak) && edgeToBreak == -1) {
    child  <- edge[, 2]
    nTips <- (length(parent) / 2L) + 1L
    samplable <- child > nTips
    # newEdges <- vapply(which(samplable), DoubleNNI, parent=parent, child=child, list(matrix(0L, nEdge, 2), matrix(0L, nEdge, 2)))
    newEdges <- unlist(lapply(which(samplable), DoubleNNI, parent=parent, child=child), recursive=FALSE) # Quicker than vapply, surprisingly
    newTrees <- lapply(newEdges, function (edges) {tree$edge <- edges; tree}) # Quicker than vapply, surprisingly
    
    # Return:
    newTrees
  } else {
    newEdge <- NNISwap(parent, edge[, 2], edgeToBreak=edgeToBreak)
    tree$edge <- cbind(newEdge[[1]], newEdge[[2]])
    
    # Return:
    tree
  }
}

#' @describeIn NNI faster version that takes and returns parent and child parameters
#' @template treeParent
#' @template treeChild
#' @param nTips (optional) Number of tips.
#' @return a list containing two elements, corresponding in turn to the rearranged parent and child parameters
#' @importFrom TreeTools SampleOne
#' @export
NNISwap <- function (parent, child, nTips = (length(parent) / 2L) + 1L, edgeToBreak=NULL) {
  rootNode  <- nTips + 1L
  samplable <- child > nTips
  if (!any(samplable)) stop("Not enough edges to allow NNI rearrangement")
  
  if (is.null(edgeToBreak)) { 
    edgeToBreak <- SampleOne(which(samplable))
  } else if (!samplable[edgeToBreak]) {
    stop("edgeToBreak must be an internal edge")
  }

  if (is.na(edgeToBreak)) stop("Cannot find a valid rearrangement")
  
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2   <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]

  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  childSwap <- child[newInd]
  child[oldInd] <- childSwap
  RenumberEdges(parent, child)
}

## TODO use RenumberList
#' Double NNI
#' 
#' Returns the edge parameter of the two trees consistent with the speficied \acronym{NNI} rearrangement
#'
#' @template treeParent
#' @template treeChild
#' @template edgeToBreakParam
#'
#' @return the \code{tree$edge} parameter of the two trees consistent with the specified rearrangement
#'
#' @keywords internal
#' @importFrom TreeTools RenumberTree
#' @author Martin R. Smith
#' 
DoubleNNI <- function (parent, child, edgeToBreak) {
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2.3 <- which(parent == end2)
  ind2   <- ind2.3[1]
  ind3   <- ind2.3[2]

  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  child2 <- child
  childSwap <- child[newInd]
  child2[oldInd] <- childSwap
  
  newInd <- c(ind3, ind1)
  oldInd <- c(ind1, ind3)
  childSwap <- child[newInd]
  child[oldInd] <- childSwap
  
  nEdge <- length(parent)
  
  # Return:
  list(RenumberTree(parent, child), RenumberTree(parent, child2))
}

#' Rooted NNI 
#' @describeIn NNI Perform \acronym{NNI} rearrangement, retaining position of root
#' @export
RootedNNI <- function (tree, edgeToBreak=NULL) {
  edge <- tree$edge
  if (!is.null(edgeToBreak) && edgeToBreak == -1) {
    parent <- edge[, 1]
    child  <- edge[, 2]
    nTips <- (length(parent) / 2L) + 1L
    rootNode <- nTips + 1L
    samplable <- parent != rootNode & child > nTips
    newEdges <- unlist(lapply(which(samplable), DoubleNNI, parent=parent, child=child), recursive=FALSE) # Quicker than vapply, surprisingly
    newTrees <- lapply(newEdges, function (edges) {tree$edge <- edges; tree}) # Quicker than vapply, surprisingly
    
    # Return:
    newTrees
  } else {
    newEdge <- RootedNNISwap(edge[, 1], edge[, 2], edgeToBreak=edgeToBreak)
    tree$edge <- cbind(newEdge[[1]], newEdge[[2]])
    
    # Return:
    tree
  }
}

#' @describeIn NNI faster version that takes and returns parent and child parameters
#' @return a list containing two elements, corresponding in turn to the rearranged parent and child parameters
#' @export
RootedNNISwap <- function (parent, child, nTips = (length(parent) / 2L) + 1L, edgeToBreak = NULL) {
  rootNode <- nTips + 1L
  
  samplable <- parent != rootNode & child > nTips

  if (is.null(edgeToBreak)) { 
    edgeToBreak <- SampleOne(which(samplable))
  } else if (!samplable[edgeToBreak]) {
    stop("edgeToBreak cannot include a tip or the root node")
  }
  
  if (is.na(edgeToBreak)) stop("Cannot find a valid rearrangement")
  
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2   <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]
  
  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  
  child_swap <- child[newInd]
  child[oldInd] <- child_swap
  RenumberEdges(parent, child)
}
