\name{ZIP}

\alias{ZIP}

\title{A zero-inflated Poisson family}

\description{
This is basically a carbon-copy of \code{gamlss.dist::ZIP()}, which establishes the zero-inflated Poisson distribution as a family for use in the main \code{rpql} function. The distribution is parametrized in terms of the mean of the Poisson distribution, and the probability of zero inflation. 

Full credit goes to the GAMLSS team for the original\code{gamlss.dist::ZIP()} family class, which was copied into this package due to issues with the \code{gamlss.dist} package being archived at some point (but might have since been resolved?!) 
}

\usage{ 
ZIP()
  }

\details{
Used in the form \code{rpql(y, ..., family = ZIP(), ...)}.
}

\value{
An object of class "family"
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}

