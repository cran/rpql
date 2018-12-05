\name{gendat.glmm}

\alias{gendat.glmm}

\title{Simulates datasets based on a Generalized Linear Mixed Model (GLMM).}

\description{
Datasets are simulated from a GLMM given a series of inputs including: model matrices \code{X} and \code{Z} for the fixed and random effects respectively, a set of true fixed effect coefficients \code{beta}, a list of true random effect covariance matrices \code{D}, the family of response, and some other nusiance parameters if appropriate.
}

\usage{ 
gendat.glmm(id, X, beta, Z, D, trial.size = 1, family = gaussian(), 
  phi = NULL, shape = NULL, zeroprob = NULL, upper.count = Inf)
   }

\arguments{
  \item{id}{A list with each element being a vector of IDs that reference the model matrix in the corresponding element in the list \code{Z}. Each vector of IDs \emph{must} be integers (but not factors).}

   \item{X}{A model matrix of corresponding to the fixed effects. A column of ones should be included if a fixed intercept is to be included in the model.}
   
  \item{beta}{A vector of true fixed effect parameters, with the same length as the number of columns in \code{X}.}

  \item{Z}{A list with each element being a model matrix for a set of random effects. Each element of \code{Z} is referenced by a vector of IDs given by the corresponding element in the list \code{id}. Each model matrix (element of \code{Z}) should have the same number of rows as the length of \code{y}.}

   \item{D}{A list with each element being a symmetric random effects covariance matrix which is used to generate random effects. These random effects are then applied to the corresponding element in the list \code{Z}, and are referenced by the corresponding element in the list \code{id}.}

  \item{trial.size}{The trial size if \code{family = binomial()}. Either takes a single non-zero value or a vector of non-zero values with length the same as the number of rows in \code{X}. The latter allows for differing trial sizes across responses. Defaults to 1.} 

  \item{family}{The distribution for the responses in GLMM. The argument must be applied as a object of class "family". Currently supported arguments include: \code{gaussian()}, \code{poisson()}, \code{binomial()}, \code{Gamma()}, \code{nb2()} for negative binomial, \code{LOGNO()} for log-normal, and \code{ZIP()} for zero-inflated Poisson.}

   \item{phi}{A non-zero value for the true variance parameter \eqn{\sigma^2} if \code{family = gaussian()}, the true variance parameter \eqn{\sigma^2} on the log scale if \code{family = LOGNO()}, or the overdispersion parameter if \code{family = nb2()}, where the negative binomial variance is parameterized as \eqn{V = \mu + \phi\mu^2}. Defaults to NULL.}

   \item{shape}{A non-zero value for the shape parameter \eqn{a} if \code{family = Gamma()}, where the variance is parameterized as \eqn{V = \mu^2/a}. Defaults to NULL.}

   \item{zeroprob}{A value between 0 and 1 for the probability of a structural zero if \code{family = ZIP()} for zero-inflated Poisson. Defaults to NULL.}
   
  \item{upper.count}{A non-zero integer which allows the user to control the maximum value of the counts generates for datasets when \code{family = poisson()} or \code{nb2()}. When the responses are simulated, a \code{while} loop is run to ensure that all responses generated are less than or equal to \code{upper.count}. Default to \code{Inf}.}
   }
   
\details{
The relationship between the mean of the responses and covariates in a GLMM is given as follows: For \eqn{i = 1,\ldots,n}, where \eqn{n} is, equivalently, the number of rows in \code{X}, the length of each element in \code{id}, and the number of rows in each element of \code{Z}, we have

\deqn{g(\mu_{i}) = \bm{x}^T_i \bm{\beta} + \bm{z}^T_{i1} \bm{b}_{i1} + \bm{z}^T_{i2} \bm{b}_{i2} + \ldots,}

where \eqn{g(\cdot)} is the link function, \eqn{\mu_i} is the mean of the distribution for observation \eqn{i}, \eqn{\bm{x}_i} is row \eqn{i} of the fixed effects model matrix \code{X}, and \eqn{\bm{\beta}} is the fixed effects coefficients. For the random effects, \eqn{\bm{z}_{i1}} is row \eqn{i} of the random effects model matrix in the first element of \code{Z}, while \eqn{\bm{b}_{i1}} is the vector of random effects generated for observation \eqn{i} based on the first element of \code{D}. The remaining parameters \eqn{\bm{z}_{i2}}, \eqn{\bm{b}_{i2}} and so on, are defined similarly.

Having lists for \code{id, Z}, and \code{D} allows for multiple sets of random effects to be included in the true GLMM. This is analogous to the \code{lme4} package, where multiple random effects are permitted in the formula, e.g., \code{(1|creek) + (1|creek:sample)}. If the true GLMM contains only one set of random effects, e.g., in longitudinal data, then the three lists will all contain only one element. Cases with multiple sets of random effects include nested and crossed designs, in which case \code{id, Z}, and \code{D} will have two or more elements. 

It is recommended that the user think through and design these lists carefully to ensure that they are actually constructing a true GLMM that they want to simulated data from. Yes it takes some getting use too, and we apologize for this =( Please see examples below for some ideas.

Finally, note that some of the elements of \code{beta} can be zero, i.e. truly unimportant fixed effects. Likewise, each element of \code{D} can be a random effects covariance matrix containing zero rows and columns, i.e. truly unimportant random effects.
}

\value{
A list containing the following elements
\item{y}{The vector simulated responses.}
\item{b}{A list with each element being a matrix of random effects simulated from a multivariate normal distribution with mean zero and covariance matrix equal to the corresponding element in the list \code{D}. For each element in \code{b}, the number of columns of the matrix equals the dimension of corresponding covariance matrix element in \code{D}, while the number of rows equals to the number of unique IDs in the corresponding element of the list \code{id}.}
\item{id, X, Z, beta, D, phi, shape, zeroprob, trial.size, family}{Some of the arguments entered into \code{gendat.glmm}.}
\item{nonzero.beta}{A vector indexing the non-zero values of \code{beta}, i.e. the truly important fixed effects.}
\item{nonzero.b}{A list with each element being a vector indexing the non-zero diagonal variances in the corresponding element of the list \code{D}, i.e. the truly important random effects.}
}

\references{
\itemize{
  \item Schielzeth, H., & Nakagawa, S. (2013). Nested by design: model fitting and interpretation in a mixed model era. Methods in Ecology and Evolution, 4, 14-24.
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

##################
## Example 1: Linear Mixed Models 
## Independent cluster model with 50 clusters
## Nine covariates including a fixed and random intercept
library(mvtnorm)
library(lme4)

n <- 50; m <- 10; p <- 8; 
## Generate rows of a model matrix from a multivariate normal distribution with 
## AR1 covariance structure. 

H <- abs(outer(1:p, 1:p, "-")) 
X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 

Z <- X 
true_betas <- c(1,3,2,1.5,-1,0,0,0,0) ## 5 important fixed effects 
true_D <- matrix(0,p+1,p+1) ## 3 important random effects
true_D[1:3,1:3] <- matrix(c(9,4.8,0.6,4.8,4,1,0.6,1,1),3,3,byrow=TRUE)

simy <- gendat.glmm(id = list(cluster = rep(1:n,each=m)), X = X, beta = true_betas, 
	Z = list(cluster = Z), D = list(cluster = true_D), phi = 1, family = gaussian()) 
## Notice how id, Z, and D all are lists with one element, and that 
## the name of the first element (a generic name "cluster") is the 
## same for all three lists. 
## id is where the action takes place. In particular, id$cluster is 
## designed so that the first 12 elements correspond to cluster 1, 
## the second 12 elements correspond to cluster 2, and so forth. 
## In turn, the first 12 rows of X and Z$cluster correspond 
## to cluster 1, and so on. 

\dontrun{
dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
fit_satlme4 <- lmer(y ~ X - 1 + (Z - 1 | cluster), data = dat,
	REML = FALSE)
fit_sat <- build.start.fit(fit_satlme4, gamma = 2)


lambda_seq <- lseq(1e-4,1,length=100)
fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
	family = gaussian(), lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)

summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
}


##################
## Example 2: Bernoulli GLMMs on simulated data
## Nested data with 200 observations in total: split into 10 creeks, 
## 5 samples nested within each creek

mn <- 200; 
X <- as.matrix(rep(1,mn)); 
ids <- list(samples = rep(1:50,each=4), creek = rep(1:10,each=20)) 
## We have two sets of random intercepts only, one for creek and one 
## 	for samples nested within creek.
Zs <- list(samples = X, creek = X) 

true_betas <- -0.1 
## Please ensure each element of true_D is a matrix
true_D <- list(samples = as.matrix(0.001), creek = as.matrix(1)) 

simy <- gendat.glmm(id = ids, X = X, beta = true_betas, Z = Zs, D = true_D, 
	trial.size = 1, family = binomial())

\dontrun{

## Construct a solution path use adaptive LASSO for selection
## Here is another way of constructing the adaptive weights:
## Use the fact that rpql can do a final fit based on maximum likelihood
## to obtain a good saturated fit.
fit_sat <- rpql(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
	family = binomial(), lambda = 0, hybrid = TRUE)
fit_sat <- build.start.fit(fit_sat$hybrid, gamma = 2)
	
lambda_seq <- lseq(1e-6,1,length=100)
fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)

summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
}
}
