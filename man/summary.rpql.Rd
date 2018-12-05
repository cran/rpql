\name{summary.rpql}

\alias{summary.rpql}

\alias{print.summary.rpql}

\title{Summary of GLMM fitted using regularized PQL.}

\description{A summary of the results from applying \code{rpql}.}

\usage{
\method{summary}{rpql}(object, ...)

\method{print}{summary.rpql}(x,...)
}

\arguments{
  \item{object}{An object of class "rpql".}
  \item{x}{An object of class "rpql".}
  \item{...}{Not used.}
  }

\value{
A list (some of which is printed) containing the following elements:
\item{Call}{The matched call.}

\item{fixed}{Estimated fixed effects coefficients.}

\item{ranef}{A list with each element being a matrix of estimated random effects coefficients.}

\item{ran.cov}{A list with each element being a estimated random effects covariance matrix.}

\item{logLik}{PQL log-likelihood value at convergence.}

\item{family}{The \code{family} argument, i.e. response type.}

\item{pen.type,lambda}{Penalties used for selection and the corresponding tuning parameter values.}

\item{ics}{A vector containing the number of estimated, non-zero parameters, and three information criterion. Please see the help file for \code{rpql} for details on these criteria.}

\item{id}{The \code{id} argument, i.e. list of IDs.}

\item{nonzero.fixef}{A vector indexing which of the estimated fixed effect coefficients are non-zero.}

\item{nonzero.ranef}{A list with each element being a vector indexing which of the estimated random effects are non-zero, i.e. which of the diagonal elements in the corresponding element of \code{ran.cov} are non-zero.}
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}

\seealso{\code{\link{rpql}} for fitting and performing model selection in GLMMs using regularized PQL.}

\examples{
## Please see examples in help file for the rpql function
}
