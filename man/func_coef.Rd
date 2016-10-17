\name{func_coef}
\alias{func_coef}

\title{
Haar wavelet coefficients
}
\description{
Computes Haar wavelet coefficients at a given scale.
}
\usage{
func_coef(z, scale)
}

\arguments{
  \item{z}{
  input time series
  }
  \item{scale}{
  a negative integer representing the Haar wavelet scale
  }
}
\details{
By padding the input series around the boundary, the returned coefficient series are of the same length as the input time series.
}

\value{
A series of Haar wavelet coefficients at a given scale.
}

\author{
Haeran Cho
}

\keyword{internal}
