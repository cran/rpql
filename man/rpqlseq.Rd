\name{rpqlseq}

\alias{rpqlseq}

\title{Wrapper function for joint effects selection in GLMMs using regularized PQL.}

\description{
\code{rpql} offers fast joint selection of fixed and random effects in Generalized Linear Mixed Model (GLMMs) via regularization. The penalized quasi-likelihood (PQL) is used as a loss function, and penalties are added on to perform fixed and random effects selection. This method of joint selection in GLMMs, referred to regularized PQL, is fast compared to information criterion and hypothesis testing (Hui et al., 2016). 

\code{rpqlseq} is a wrapper function to permit a sequence of tuning parameter values, which wraps around the code workhorse function \code{rpql}.
}

\usage{
rpqlseq(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, 
  pen.type = "lasso", start = NULL, cov.groups = NULL, pen.weights = NULL, 
  offset = NULL, intercept = TRUE, save.data = FALSE, 
  control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, 
  scad.a = 3.7, mcp.gamma = 2, seed = NULL), ...)

  }

  
\arguments{
  \item{y, X, Z, id, family, trial.size}{As per the \code{rpql} function. Please see the help file for \code{rpql} for details on the arguments.}

  \item{lambda}{Either a vector containing \bold{sequence} of tuning parameter values, which is applied to both penalties, or two-column matrix containing a \bold{sequence} of tuning parameter values for the fixed and random effects penalty respectively.}

  \item{pen.type, start, cov.groups, pen.weights, offset, intercept, save.data, control}{As per the \code{rpql} function. Please see the help file for \code{rpql} for details on the arguments.}

  \item{...}{Not used.}
}

\details{
Please see the help file for \code{rpql} for details on how regularized PQL works. \code{rpqlseq} is simply a wrapper function to run the core \code{rpql} function multiple times, on a sequence of tuning parameter values, in order to construct a regularization path. The best models, based on different information criteria for selecting the best tuning parameter (degree of sparsity) are then returned.
}

\value{
An object of class "rpql" containing the following elements:
\item{best.fits}{A list containing the best fitted models as based on different information criteria used to select the tuning parameter. Each element in this list has the same structure as the output from the \code{rpql} function. Please see the \code{rpql} function for details on the information criteria available as well as the nature of the output.}

\item{collect.ics}{A matrix containing the values of various information criteria calculated for the sequence of \code{lambda} values supplied. The best fitted models found in \code{best.fits} is based off this matrix i.e., each element in \code{best.fits} corresponds to a model that was chosen based on minimizing the corresponding information criterion in \code{collect.ics}. Please see the \code{rpql} function for details on the information criteria available.}

\item{lambda}{The sequence of tuning parameters considered.}
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}


\seealso{
\code{\link{rpql}}, which is the core workhorse function that performed regularized PQL for a single set of tuning parameter values.
}

\examples{
## Please see examples in help file for the \code{rpql} function for usage.
}
