% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SuccessiveApproximations.R
\name{Suboptimality}
\alias{Suboptimality}
\title{Tree suboptimality}
\usage{
Suboptimality(trees, proportional = FALSE)
}
\arguments{
\item{trees}{list of trees, to include an optimal tree}

\item{proprtional}{logical stating whether to normalise results to lowest score}
}
\value{
a vector listing, for each tree, how much their score differs from the optimal (lowest) score.
}
\description{
How suboptimal is a tree?
}
\keyword{internal}
