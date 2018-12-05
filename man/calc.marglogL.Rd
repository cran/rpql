\name{calc.marglogL}

\alias{calc.marglogL}

\title{Calculate the marginal log-likelihood for a GLMM fitted using \code{rpql}}

\description{
After fitting and performing joint (fixed and random effects) using regularized PQL, one may then (for one reason or another) want to calculate the marginal likelihood for the (sub)model, possibly on a test dataset for prediction. This is the main purpose of \code{calc.marglogL}.
}

\usage{
calc.marglogL(new.data, fit, B = 1000)
  }

\arguments{
  \item{new.data}{A list containing the elements \code{new.data$y},  \code{new.data$X},  and \code{new.data$Z}. These correspond respectively to the responses, fixed effects model matrix, and random effects model matrices that the marginal log-likelihood is be calculated on. No check is made against the elements in \code{fit} to ensure that these are of the correct dimensions compared, and furthermore it is assumed that \code{new.data$Z} is a list in the same order as the \code{Z} used when fitting the original model via \code{rpql}.}

  \item{fit}{An object of class \code{pqrl}. In the least, \code{fit} should be a list containing the elements \code{fit$family} for the family, e.g. gaussian(), poisson(), \code{fit$fixef} for the estimated vector of fixed effects, \code{fit$ran.cov} which is a list of estimated random effects covariance matrices. If appropriate, \code{fit} may also contain the elements \code{fit$phi} for the estimated variance parameter in normal, lognormal, and negative binomial GLMMs, \code{fit$shape} for the estimated shape parameter used in Gamma GLMMs, \code{fit$trial.size} for the trial size(s) for binomial GLMMs, and \code{fit$zeroprob} for the estimated probability of a structural zero in ZIP GLMMs.}

  \item{B}{A positive integer for the number of random effects examples to generate, when performing Monte-Carlo integration. Defaults to 1000.}
  }
  
  
\details{
Regularized PQL performs penalized joint (fixed and random effects) selection for GLMMs, where the penalized quasi-likelihood (PQL, Breslow and Clayton, 1993) is used the loss function. After fitting, one may then wish to calculate the marginal log-likelihood for the (sub)model, defined as 

\deqn{\ell = \log\left(\int f(\bm{y}; \bm{\beta}, \bm{b}, \phi) f(\bm{b}; \bm{\Sigma}) d\bm{b}\right),}

where \eqn{f(\bm{y}; \bm{\beta}, \bm{b}, \phi)} denotes the conditional likelihood of the responses \eqn{\bm{y}} given the fixed effects \eqn{\bm{\beta}}, random effects \eqn{\bm{b}}, and nuisance parameters \eqn{\phi} if appropriate, and \eqn{f(\bm{b}; \bm{\Sigma})} is the multivariate normal distribution for the random effects, with covariance matrix \eqn{\bm{\Sigma}}. \code{calc.marglogL} calculates the above marginal likelihood using Monte-Carlo integration.

Admittedly, this function is not really useful for fitting the GLMM \emph{per-se}: it is never called by the main function \code{rpql}, and the marginal likelihood is (approximately) calculated anyway if \code{hybrid.est = TRUE} and the final submodel is refitted using \code{lme4}. Where the function comes in handy is if you have a validation or test dataset, and you want to calculated the predicted (log) likelihood of the test data given the regularized PQL fit.  
}

\value{The marginal log-likelihood of \code{new.data} given the GLMM in \code{fit}.}

\references{
\itemize{
\item Breslow, N. E., & Clayton, D. G. (1993). Approximate inference in generalized linear mixed models. Journal of the American Statistical Association, 88, 9-25.
  }
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}


\section{Warnings}{
  \itemize{
  \item No check is made to see if the dimensions of the elements \code{new.data} and \code{fit} match, e.g. the number of columns in \code{new.data$X} is equal to the number of elements in \code{fit$fixef}. Please ensure they are!
  \item Monte-Carlo integration is computationally intensive especially if \eqn{\bm{y}} is long!
    } 
  }


\seealso{
\code{\link{rpql}} for fitting and performing model selection in GLMMs using regularized PQL. \code{lme4} also approximately calculates the marginal log-likelihood when fitting a GLMM.}

\examples{
## Make your own =D
}
