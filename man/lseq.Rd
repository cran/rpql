\name{lseq}

\alias{lseq}

\title{Generates a sequence of tuning parameters on the log scale}

\description{
Generates a sequence of tuning parameters \eqn{\lambda} that are equally spaced on the log-scale. It may be used as part of constructing a solution path for the main fitting function \code{rpql}.
}

\usage{
lseq(from, to, length, decreasing = FALSE)
  }

\arguments{
  \item{from}{The minimum tuning parameter to start the sequence from.}
  \item{to}{The maximum tuning parameter to go to.}
  \item{length}{The length of the sequence.}
  \item{decreasing}{Should the sequence be in ascending or descending order?}
  }

\details{
  For joint selection of fixed and random effects in GLMMs, regularized PQL (Hui et al., 2016) works taking the penalized quasi-likelihood (PQL, Breslow and Clayton, 1993) as a loss function, and then sticking on some penalties in order to model variable. The penalties will depend upon one or more tuning parameters \eqn{\lambda > 0}, and the typical way this is chosen is to construct a sequence of \eqn{\lambda} values, fit the regularized PQL to each one value, and then use a method like information criterion to select the best \eqn{\lambda} and hence the best model. Please see the help file for \code{\link{rpql}} for more details, and \code{glmnet} (Friedman et al., 2010) and \code{ncvreg} (Breheny, and Huang, 2011) as examples of other packages that do penalized regression and involve tuning parameter selection.

The idea of equally spacing the sequence of \eqn{\lambda}'s on the log (base 10) scale may not necessary be what you want to do, and one is free to use the standard \code{seq()} function for constructing sequences. By equaling spacing them on log-scale, it means that there will be a large concentration of small tuning parameter values, with less large tuning parameter values (analogous to a right skewed distribution). This may be useful if you believe the that most of the penalization/variable selection action takes place on smaller values of \eqn{\lambda}.

It is somewhat of an art form to construct a good sequence of tuning parameter values: the smallest \eqn{\lambda} should produce the saturated model if possible, and the largest \eqn{\lambda} should shrink most if not all covariates to zero i.e., the null model. Good luck!
    }

\value{A sequence of tuning parameter values of length equal to \code{length}.}

\references{
\itemize{
\item Breheny, P. and Huang, J. (2011) Coordinate descent algorithms fof nonconvex penalized regression, with applications to biological feature selection. The Annals of Appliedv Statistics, 5, 232-253.
\item Breslow, N. E., and Clayton, D. G. (1993). Approximate inference in generalized linear mixed models. Journal of the American Statistical Association, 88, 9-25.
\item Friedman, J., Hastie T., and Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33, 1-22. URL: http://www.jstatsoft.org/v33/i01/.
   \item Hui, F.K.C., Mueller, S., and Welsh, A.H. (2016). Joint Selection in Mixed Models using Regularized PQL. Journal of the American Statistical Association: accepted for publication.
  }
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}


\seealso{
\code{\link{rpql}} for fitting and performing model selection in GLMMs using regularized PQL.
}


\examples{
## Please see examples in help file for the rpql function
}
 

  
  
 