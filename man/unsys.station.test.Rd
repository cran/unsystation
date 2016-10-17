\name{unsys.station.test}
\alias{unsys.station.test}

\title{
A second-order stationarity of time series based on unsystematic sub-samples
}

\description{
The function implements a stationarity test procedure, where the main statistic is obtained from measuring the difference in the second-order structure over pairs of randomly drawn intervals. Maximising the main statistics after AR Sieve bootstrap-based stabilisation, the test statistic is drawn, which is reported along with the corresponding pair of intervals and the test result.
}

\usage{
unsys.station.test(x, M = 2000, sig.lev = 0.05, max.scale = NULL, m =
                 NULL, B = 200, eps = 5, use.all = FALSE, do.parallel =
                 FALSE)
}

\arguments{
  \item{x}{
  input time series
  }
  \item{M}{
  number of randomly drawn intervals
  }
  \item{sig.lev}{
  significance level between \code{0} and \code{1}
  }
  \item{max.scale}{
  number of wavelet scales used for wavelet periodogram computation
  }
  \item{m}{
  minimum length of a random interval
  }
  \item{B}{
  bootstrap sample size
  }
  \item{eps}{
  a parameter used for random interval generation, see the supplementary document of Cho (2016)
  }
  \item{use.all}{
  if \code{use.all=TRUE}, all \code{M*M} pairs of random intervals are considered in test statistic computation; 
  if \code{use.all=FALSE}, only \code{10*M} pairs are used; regardless, the whole \code{M*M} pairs are considered 
  in test criterion generation
  }
  \item{do.parallel}{
  if \code{do.parallel=TRUE}, a set of copies of R running in parallel are created and used for bootstrap procedure
  }
}

\value{
  \item{intervals}{
  a pair of intervals corresponding to the test statistic, exhibiting the most distinct second-order behaviour
  }
  \item{test.stat}{
  test statistic
  }
  \item{test.criterion}{
  test criterion
  }
  \item{test.res}{
  if \code{test.res=TRUE}, the null hypothesis of stationarity is rejected at the given significance level
  }
}

\references{
H. Cho (2016) A second-order stationarity of time series based on unsystematic sub-samples. Stat (to appear).
}

\author{
Haeran Cho
}

\examples{
x <- rnorm(200)
unsys.station.test(x, M=1000)
}
