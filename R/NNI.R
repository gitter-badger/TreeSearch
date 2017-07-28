#' @title Tree rearrangement functions
#' 
#' These functions performs a single random \acronym{TBR}, \acronym{SPR} or \acronym{NNI} iteration.
#'
#' Performs a single iteration of the nearest-neigbour interchange, subtree pruning and regrafting,
#' or tree bisection and reconnection algorithms.
#' NNI and SPR are based on the corresponding phangorn functions, but have been re-coded to 
#' improve their speed.
#' 
#' Branch lengths are not supported.
#' 
#' @return Returns a tree with class \code{phylo}.
#'
#' @template treeParam
#' @param edgeToBreak the index of an edge to bisect, generated randomly if not specified
#' 
#' @references
#' The algorithms are summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' @author Martin R. Smith
#' 
#' @examples
#' tree <- ape:::rtree(20, br=NULL)
#' NNI(tree)
#' SPR(tree)
#' TBR(tree)
#' @export
NNI <- function (tree) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nTips  <- length(tree$tip.label)
  rootNode <- nTips + 1L
  chosenInternalEdge <- SampleOne(which(child > nTips))
  if(is.na(chosenInternalEdge)) return(NULL)
  nEdge <- length(parent)
  nNode <- tree$Nnode
  if (nNode == 1) return(tree)
  
  end1  <- parent[chosenInternalEdge]
  end2  <- child[chosenInternalEdge]
  ind1    <- which(parent == end1)
  ind1    <- ind1[ind1 != chosenInternalEdge][1]
  ind2    <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]
  new_ind <- c(ind2, ind1)
  old_ind <- c(ind1, ind2)
  child_swap <- child[new_ind]
  child[old_ind] <- child_swap
  tree$edge <- OrderEdgesNumberNodes(parent, child, nTips, nEdge)
  tree
}

#' Rearrange a rooted tree
#'
#' This function performs a rearrangement iteration on a tree, retaining the position of the root.
#'
#' A single \acronym{NNI}, \acronym{SPR} or \acronym{TBR} rearrangement is performed, subject to the constraint that 
#' no taxon may be moved to the opposite side of the root node.
#' Branch lengths are not (yet) supported.
#' 
#' @usage
#' RootedNNI(tree)
#' RootedSPR(tree)
#' RootedTBR(tree)
#'
#' @param tree A bifurcating tree of class \code{\link{phylo}}, with all nodes resolved;
#' 
#' @return This function returns a tree, in \code{phylo} format.
#'
#' @author Martin R. Smith
#' \code{RootedNNI} is abridged from the \pkg{phangorn} function \code{nnin}
#' 
#' @seealso
#' \itemize{
#' \item \code{\link{SetOutgroup}}, set the outgroup of the phylogenetic tree
#' \item \code{\link{NNI}}, unrooted \acronym{NNI} and \acronym{SPR}
#' \item \code{\link{TBR}}, unrooted \acronym{TBR}
#' }
#' 
#' @examples{
#'   require('ape')
#'   tree <- read.tree(text='(((a,b),c),(d,(e,f)));')
#'   tree <- SetOutgroup(tree, c('e', 'f'))
#'   plot(tree)
#'   dev.new()
#'   plot(RootedNNI(tree))
#'   plot(RootedSPR(tree))
#'   plot(RootedTBR(tree))
#' }
#' 
#'
#' @export
RootedNNI <- function (tree) {
  nTips  <- length(tree$tip.label)
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  sampleableChild <- child
  sampleableChild[which(parent == as.integer(parent[!match(parent, child, 0)][1]))] <- -1 # Don't want to switch across the root
  chosenInternalEdge <- SampleOne(which(sampleableChild > nTips))
  if(is.na(chosenInternalEdge)) return(NULL)
  
  rootNode <- nTips + 1L
  nEdge <- length(parent)
  nNode <- tree$Nnode
  if (nNode == 1) return(tree)
  
  end1  <- parent[chosenInternalEdge]
  end2  <- child[chosenInternalEdge]
  ind1    <- which(parent == end1)
  ind1    <- ind1[ind1 != chosenInternalEdge][1]
  ind2    <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]
  new_ind <- c(ind2, ind1)
  old_ind <- c(ind1, ind2)
  child_swap <- child[new_ind]
  child[old_ind] <- child_swap
  tree$edge <- OrderEdgesNumberNodes(parent, child, nTips, nEdge)
  tree
}