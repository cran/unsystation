\name{funcRes}
\alias{funcRes}

\title{
Random interval generation and main statistic computation
}

\description{
Called for within \code{\link{unsys.station.test}}, the function generates random intervals and computes the main statistics over pairs of the randomly drawn intervals.
}

\usage{
funcRes(yMat, M, minLength, probVec, topCand, bottomCand, var)
}

\arguments{
  \item{yMat}{
  a matrix containing wavelet periodograms at multiple scales as its columns
  }
  \item{M}{
  number of randomly drawn intervals
  }
  \item{minLength}{
  minimum length of a random interval
  }
  \item{probVec}{
  a vector assigning equal probabilities for all elements in \code{topCand} and \code{bottomCand}
  }
  \item{topCand, bottomCand}{
  candidate `positions' for randomly drawn intervals, see the supplementary document of Cho (2016).
  }
  \item{var}{
  a vector containing scale-wise average of wavelet periodograms of the input time series for standardisation
  }
}

\value{
  \item{res}{
  a matrix containing the start and end points of the pairs of random intervals, 
  and the corresponding main statistics at multiple scales
  }
  \item{se}{
  a matrix containing the start and end points of the pairs of random intervals
  }
}

\references{
H. Cho (2016) A second-order stationarity of time series based on unsystematic sub-samples. Stat (to appear).
}

\author{
Haeran Cho
}

\keyword{internal}
