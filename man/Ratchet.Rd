% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Jackknife.R, R/Ratchet.R
\name{Jackknife}
\alias{Jackknife}
\alias{Ratchet}
\alias{ProfileRatchet}
\alias{IWRatchet}
\alias{MultiRatchet}
\alias{IWMultiRatchet}
\alias{RatchetConsensus}
\alias{IWRatchetConsensus}
\title{Parsimony Ratchet}
\usage{
Jackknife(tree, dataset, resampleFreq = 2/3,
  InitializeData = PhyDat2Morphy, CleanUpData = UnloadMorphy,
  TreeScorer = MorphyLength, EdgeSwapper = TBRSwap, jackIter = 5000L,
  searchIter = 4000L, searchHits = 42L, verbosity = 1L, ...)

Ratchet(tree, dataset, InitializeData = PhyDat2Morphy,
  CleanUpData = UnloadMorphy, TreeScorer = MorphyLength,
  Bootstrapper = MorphyBootstrap, swappers = list(TBRSwap, SPRSwap,
  NNISwap), BootstrapSwapper = if (class(swappers) == "list")
  swappers[[length(swappers)]] else swappers, returnAll = FALSE,
  stopAtScore = NULL, stopAtPeak = FALSE, stopAtPlateau = 0L,
  ratchIter = 100, ratchHits = 10, searchIter = 4000,
  searchHits = 42, bootstrapIter = searchIter,
  bootstrapHits = searchHits, verbosity = 1L, suboptimal = 1e-08,
  ...)

ProfileRatchet(tree, dataset, swappers = list(TBRSwap, SPRSwap, NNISwap),
  BootstrapSwapper = if (class(swappers) == "list")
  swappers[[length(swappers)]] else swappers, returnAll = FALSE,
  stopAtScore = NULL, stopAtPeak = FALSE, stopAtPlateau = 0L,
  ratchIter = 100, ratchHits = 10, searchIter = 2000,
  searchHits = 40, bootstrapIter = searchIter,
  bootstrapHits = searchHits, verbosity = 1L, suboptimal = 1e-08,
  ...)

IWRatchet(tree, dataset, concavity = 10, swappers = list(TBRSwap,
  SPRSwap, NNISwap), BootstrapSwapper = if (class(swappers) == "list")
  swappers[[length(swappers)]] else swappers, returnAll = FALSE,
  stopAtScore = NULL, stopAtPeak = FALSE, stopAtPlateau = 0L,
  ratchIter = 100, ratchHits = 10, searchIter = 2000,
  searchHits = 40, bootstrapIter = searchIter,
  bootstrapHits = searchHits, verbosity = 1L, suboptimal = 1e-08,
  ...)

MultiRatchet(tree, dataset, ratchHits = 10, searchIter = 500,
  searchHits = 20, verbosity = 0L, swappers = list(RootedNNISwap),
  nSearch = 10, stopAtScore = NULL, ...)

IWMultiRatchet(tree, dataset, ratchHits = 10, concavity = 4,
  searchIter = 500, searchHits = 20, verbosity = 0L,
  swappers = list(RootedNNISwap), nSearch = 10,
  suboptimal = suboptimal, stopAtScore = NULL, ...)

RatchetConsensus(tree, dataset, ratchHits = 10, searchIter = 500,
  searchHits = 20, verbosity = 0L, swappers = list(RootedNNISwap),
  nSearch = 10, stopAtScore = NULL, ...)

IWRatchetConsensus(tree, dataset, ratchHits = 10, concavity = 4,
  searchIter = 500, searchHits = 20, verbosity = 0L,
  swappers = list(RootedNNISwap), nSearch = 10,
  suboptimal = suboptimal, stopAtScore = NULL, ...)
}
\arguments{
\item{tree}{A tree of class \code{\link{phylo}}.}

\item{dataset}{a dataset in the format required by \code{TreeScorer()}.}

\item{resampleFreq}{Double between 0 and 1 stating proportion of characters to resample}

\item{InitializeData}{Function that sets up data object to prepare for tree search.
The function will be passed the \code{dataset} parameter.
Its return value will be passed to \code{TreeScorer()} and \code{CleanUpData()}.}

\item{CleanUpData}{Function to destroy data object on function exit.
The function will be passed the value returned by \code{InitializeData()}.}

\item{TreeScorer}{function to score a given tree.
The function will be passed three parameters, corresponding to the
\code{parent} and \code{child} entries of a tree's edge list, and a dataset.}

\item{EdgeSwapper}{a function that rearranges a parent and child vector,
and returns a list with modified vectors; for example \code{\link[=SPRSwap]{SPRSwap()}}.}

\item{jackIter}{Integer specifying number of jackknife iterations to conduct}

\item{searchIter}{Integer specifying maximum rearrangements to perform on each bootstrap or
ratchet iteration.
To override this value for a single swapper function, set e.g.
\code{attr(SwapperFunction, 'searchIter') <- 99}}

\item{searchHits}{Integer specifying maximum times to hit best score before terminating a tree
search within a ratchet iteration.
To override this value for a single swapper function, set e.g.
\code{attr(SwapperFunction, 'searchHits') <- 99}}

\item{verbosity}{Level of detail to display in console: larger numbers provide
more verbose feedback to the user.}

\item{\dots}{further arguments to pass to \code{TreeScorer()}
(e.g. \code{\link{TipsAreColumns}}, \code{dataset}).}

\item{Bootstrapper}{Function to perform bootstrapped rearrangements of tree.
First arguments will be an \code{edgeList} and a dataset, initialized using
\code{InitializeData()}. Should return a rearranged \code{edgeList}.}

\item{swappers}{A list of functions to use to conduct edge rearrangement during tree search.
Provide functions like \code{\link{NNISwap}} to shuffle root position,
or \code{\link{RootedTBRSwap}} if the position of the root should be retained.
You may wish to use extreme swappers (such as \acronym{TBR}) early in the list,
and a more subtle rearranger (such as \acronym{NNI}) later in the list to make
incremental tinkerings once an almost-optimal tree has been found.}

\item{BootstrapSwapper}{Function such as \code{\link{RootedNNISwap}} to use
to rearrange trees within \code{Bootstrapper()}.}

\item{returnAll}{Set to \code{TRUE} to report all MPTs encountered during the
search, perhaps to analyse consensus.}

\item{stopAtScore}{stop search as soon as this score is hit or beaten.}

\item{stopAtPeak}{Logical specifying whether to terminate search once a subsequent
iteration recovers a sub-optimal score. Useful with methods that
return all trees one rearrangement from the current tree, such
as \code{\link[=AllTBR]{AllTBR()}}.  Will be overridden if a passed function
has an attribute \code{stopAtPeak} set by
\code{attr(FunctionName, 'stopAtPeak') <- TRUE}.}

\item{stopAtPlateau}{Integer. If > 0, tree search will terminate if the score
has not improved after \code{stopAtPlateau} iterations.
Useful with methods that return all trees one rearrangement from the current
tree, such as \code{\link[=AllTBR]{AllTBR()}}.
Will be overridden if a passed function has an attribute \code{stopAtPlateau} set
by \code{attr(FunctionName, 'stopAtPlateau') <- TRUE}.}

\item{ratchIter}{Stop when this many ratchet iterations have been performed.}

\item{ratchHits}{Stop when this many ratchet iterations have found the same
best score.}

\item{bootstrapIter}{Integer specifying maximum rearrangements to perform on each bootstrap
iteration (default: \code{searchIter}).}

\item{bootstrapHits}{Integer specifying maximum times to hit best score on each bootstrap
iteration (default: \code{searchHits}).}

\item{suboptimal}{retain trees that are suboptimal by this score. Defaults to 1e-08 to counter rounding errors.}

\item{concavity}{A numeric value to use as the concavity constant (\code{k})
in implied weighting. Defaults to 10; TNT sets a default of 3, but this
is too low in some circumstances (Smith, 2019).}

\item{nSearch}{Number of Ratchet searches to conduct
(for \code{RatchetConsensus()})}
}
\value{
a list of trees recovered after jackknife iterations

This function returns a tree modified by parsimony ratchet iterations.

\code{MultiRatchet()} returns a list of optimal trees
produced by \code{nSearch}
ratchet searches, from which a consensus tree can be generated using
\code{\link[ape:consensus]{ape::consensus()}} or \code{\link[TreeTools:ConsensusWithout]{TreeTools::ConsensusWithout()}}.

\code{IWMultiRatchet} returns a list of optimal trees produced by \code{nSearch}
Ratchet searches, using implied weighting.
}
\description{
\code{Ratchet} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
}
\section{Functions}{
\itemize{
\item \code{Jackknife}: Jackknife resampling. Note that at present this assumes that
\code{InitializeData} will return a morphy object; if this doesn't hold for you, please
let me know and I'll make the function more general.

\item \code{ProfileRatchet}: Shortcut for Ratchet search under Profile Parsimony

\item \code{IWRatchet}: Shortcut for Ratchet search using implied weights

\item \code{RatchetConsensus}: deprecated alias for \code{MultiRatchet()}

\item \code{IWRatchetConsensus}: deprecated alias for \code{MultiRatchet()}
}}

\examples{
data('Lobo', package='TreeTools')
njtree <- TreeTools::NJTree(Lobo.phy)
# Increase value of ratchIter and searchHits to do a proper search
quickResult <- Ratchet(njtree, Lobo.phy, ratchIter = 2, searchHits = 3)

# Plot result (legibly)
oldPar <- par(mar = rep(0, 4), cex = 0.75)
plot(quickResult)
par(oldPar)

# IW search is currently much slower:
quickIWResult <- IWRatchet(quickResult, Lobo.phy, concavity = 2.5,
                           ratchIter = 1, searchIter = 25, searchHits = 2,
                           swappers = RootedTBRSwap, verbosity = 5)
 
}
\references{
\itemize{
\item \insertRef{Nixon1999}{TreeSearch}
\item \insertRef{SmithTern}{TreeSearch}
}
}
\seealso{
\itemize{
\item \code{\link[=JackLabels]{JackLabels()}}: Label nodes of a tree with jackknife supports.
}

\itemize{
\item \code{\link[=TreeSearch]{TreeSearch()}}
\item Adapted from \code{\link[phangorn:parsimony]{pratchet}} in the
\pkg{phangorn} package.
}
}
\author{
Martin R. Smith

Martin R. Smith
}
\keyword{tree}
