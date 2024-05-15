\name{mlelod-package}
\alias{mlelod-package}
\docType{package}
\title{
\packageTitle{mlelod}
}
\description{
\packageDescription{mlelod}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{mlelod}
\packageIndices{mlelod}
Values below the limit of detection (LOD) are a problem in several fields of science, and there are
numerous approaches for replacing the missing data. Thic package uses a new mathematical solution for
maximum likelihood estimation that allows us to estimate the true values of the mean and standard deviation
for normal distributions and is significantly faster than previous implementations. The core function
is the function mlelod with three parameters: the size of the sample, the values of the sample and the
value of limit of detection. The function returns two estimates: for mu (expected value) and sigma
(standard deviation).
}
\author{
\packageAuthor{mlelod}

Maintainer: \packageMaintainer{mlelod}
}
\references{
The article with the derived method is submited to Journal of Statistical Software
}
\keyword{ package }
