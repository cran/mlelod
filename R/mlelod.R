mlelod <-
function (n, censoreddata, lod) 
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
