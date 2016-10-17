\name{funcApplyVar}
\alias{funcApplyVar}

\title{
Computation of the standard errors of main statistics over the bootstrap samples
}

\usage{
funcApplyVar(nullStatMat, maxLevel, Msq)
}

\description{
Computes the standard errors of main statistics over the bootstrap samples.
}

\arguments{
  \item{nullStatMat}{
  a matrix stacking the main statistics over bootstrap samples
  }
  \item{maxLevel}{
  number of columns for a resulting matrix
  }
  \item{Msq}{
  number of rows for a resulting matrix
  }
}

\value{
A matrix containing standard errors of the main statistics computed at \code{maxLevel} scales over \code{Msq} pairs of intervals.
}

\author{
Haeran Cho
}

\keyword{internal}
