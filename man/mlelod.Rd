\name{mlelod}
\alias{mlelod}
\title{
Estimates the parameters of the normal distribution.
}
\description{
The function returns the two estimates for the parameters of the normal distribution.
}
\usage{
mlelod(n, censoreddata, lod)
}
\arguments{
  \item{n}{
n is the size of the sample
}
  \item{censoreddata}{
censoreddata is the vector containing the values of the sample which are greater than lod}
  \item{lod}{
lod is the value of level of detection.
}
}
\value{
\item{muEst }{Description of 'comp1'}
\item{sigmaEst }{Description of 'comp2'}
}
\references{
The article with the derived method is submited to Journal of Statistical Software
}
\author{
Gregor Sega
}
\examples{
## The function is currently defined as
mlelod <- function (n, censoreddata, lod) 
{
    k <- length(censoreddata)
    s <- sum(censoreddata)
    s2 <- sum(censoreddata^2)
    g <- function(x, n, k, s, s2, lod) {
        a <- (k * (lod^2 - x^2) + s2 - 2 * lod * s)/(x * (k * 
            lod - s))
        b <- s - k * (k * x^2 + lod * s - s2)/(k * lod - s) - 
            x * (n - k) * dnorm(a, 0, 1)/pnorm(a, 0, 1)
        return(b)
    }
    h <- function(x) {
        return(g(x, n, k, s, s2, lod))
    }
    lsi <- sqrt(s2/k - (s/k)^2)
    sigmaEst <- uniroot(h, lower = lsi/4, upper = 4 * lsi, tol = 1e-06)$root
    muEst <- (k * sigmaEst^2 + lod * s - s2)/(k * lod - s)
    return(list(muEst = muEst, sigmaEst = sigmaEst))
  }
##define the parameters of the normal distribution
mu <- 5
sigma <- 4
##define the size of the sample and the value of lod
n <- 100
lod <- 2
##generate normally distributed data and extract the observable values (the ones exceeding lod)
data <- rnorm(n, mu, sigma)
data2 <-  Filter(function(x) x>lod,data)
##run the function to obtain the estimates
mlelod(n, data2, lod)
}
\keyword{ univar }
\keyword{ models }
