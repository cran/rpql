\name{nb2}

\alias{nb2}

\title{A negative binomial family}

\description{
Since the negative binomial is not a family in base \code{R}, an \code{nb2()} family has been created which establishes the negative binomial as a family for use in the main \code{rpql} function. Only the log link is available at the moment, with the variance parameterized as \eqn{V = \mu + \phi\mu^2} where \eqn{\phi} is the overdispersion parameter.
}

\usage{ 
nb2()
  }

\details{
Used in the form \code{rpql(y, ..., family = nb2(), ...)}.
}

\value{
An object of class "family"
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}

\seealso{
\code{\link{negative.binomial}} in the \code{MASS} package for another example of a negative.binomial family. 
}

\examples{
\dontrun{
## The function is currently defined as follows
nb2 <- function () {
    link <- "log"
    linkfun <- function(mu) log(mu)
    linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
    mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
    variance <- function(mu, phi) mu + phi * mu^2
    valideta <- function(eta) TRUE
    validmu <- function(mu) all(mu > 0)
    structure(list(family = "negative.binomial", link = "log", 
        linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        variance = variance, valideta = valideta, validmu = validmu, 
        name = link), class = "family")
  }
}  
}
