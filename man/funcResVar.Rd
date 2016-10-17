\name{funcResVar}
\alias{funcResVar}

\title{
Main statistic computation for bootstrap samples
}

\description{
Called for within \code{\link{unsys.station.test}}, the function computes the main statistics for bootstrap samples over the pairs of given intervals.
}

\usage{
funcResVar(yMat, seMat, tmpVar)
}

\arguments{
  \item{yMat}{
  a matrix containing wavelet periodograms of the bootstrapped time series
  }
  \item{seMat}{
  a matrix containing the pairs of randomly drawn intervals
  }
  \item{tmpVar}{
  a vector containing scale-wise averages of wavelet periodograms of the bootstrapped time series for standardisation
  }
}

\value{
A matrix containing the main statistics computed at multiple scales over the given pairs of intervals.
}

\author{
Haeran Cho
}

\keyword{internal}
