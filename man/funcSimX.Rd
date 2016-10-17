\name{funcSimX}
\alias{funcSimX}

\title{
AR series generation
}
\description{
The function takes AR coefficients and a matrix of innovations and generates multiple AR process.
}

\usage{
funcSimX(coef, buffMat)
}

\arguments{
  \item{coef}{
  a vector of AR coefficients
  }
  \item{buffMat}{
  a matrix of innovations including `buffers'
  }
}

\author{
Haeran Cho
}

\value{
When \code{buffMat} is \code{(T+p)}-by-\code{B} and \code{coef} is of length \code{p},
\code{B} AR series of length \code{T} are returned as the columns.
}

\examples{
coef <- .5
boot.x <- funcSimX(coef, matrix(rnorm(501*200), ncol=200))
}

\keyword{internal}
