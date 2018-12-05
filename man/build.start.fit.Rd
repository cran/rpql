\name{build.start.fit}

\alias{build.start.fit}

\title{Constructs a start fit for use in the \code{rpql} function}

\description{Takes a GLMM fitted using the \code{lme4} package i.e., using either the \code{lmer} or \code{glmer} functions, and construct a list containing starting values for use in the \code{start} argument in main fitting function \code{rpql}. It also constructs adaptive lasso weights, which can subsequently be used in the \code{pen.weights} arguments in the \code{rpql} function, if the adaptive lasso penalty is used for variable selection.}

\usage{
build.start.fit(lme4.fit, id = NULL, gamma = 0, cov.groups = NULL)
}

\arguments{
  \item{lme4.fit}{An object of class "lmerMod" or "glmerMod", obtained when fitting a (G)LMM using the \code{lmer} and \code{glmer} functions in the \code{lme4} package.}
  
  \item{id}{A optional list with each element being a vector of IDs that reference the model matrix in the corresponding element in the list \code{Z}. Each vector of IDs \emph{must} be integers (but not factors). Note this is optional argument as it is only use for non-compulsory formatting purposes in the function.}

  \item{gamma}{A vector of power parameters, \eqn{\gamma}, for use in constructing adaptive lasso weights. Can be a vector of one or two elements. If two elements, then the first and second elements are the power parameter for the fixed and random effect weights respectively. If one element, the same power parameter is used for both fixed and random effect weights. Defaults to 0, in which case the weights are all equal to 1 i.e., it reverts to the unweighted lasso penalty.}

  \item{cov.groups}{A vector specifying if fixed effect coefficients (including the intercept) should be regarded and therefore penalized in groups. For example, if one or more of the fixed effect covariates are factors, then \code{lme4} will automatically create dummy variables in the model matrix and estimate coefficients for each level, using one level as the reference. \code{cov.groups} is then used to identify all the coefficients that corresponds to that factor, such that all of these coefficients are penalized collectively as a group. Defaults to NULL, in which case it is assumed all coefficients should be treated independently. Please see the details and examples for more details.}
  }
  
\details{
This function is mainly used when: 1) you want to produce good starting values for the main fitting function \code{rpql}, and so you fit a saturated (full) GLMM using \code{lme4} and use the estimates from there as starting values, and/or 2) you want to obtain adaptive lasso weights of the form \eqn{weight_k = |\tilde{parameter}_k|^{-\gamma}}, where \eqn{\gamma > 0} is the power parameter and \eqn{\tilde{parameter}_k} is the parameter estimate from the saturated model fit. For regularized PQL specifically, this function will construct adaptive lasso weights from the \code{lme4} fit as follows: Let \eqn{w^F} and \eqn{w^R} denote fixed and random effect adaptive weights respectively. Then we have,

\deqn{w^F_k = |\tilde{\beta}_k|^{-\gamma_1}}
\deqn{w^R_l = |\tilde{\Sigma}_{ll}|^{-\gamma_2},}

where \eqn{\tilde{\beta}_k} is the estimated coefficient for the \eqn{k^{th}} fixed effect, \eqn{\tilde{\Sigma}_{ll}} is the \eqn{l^{th}} diagonal element from the estimated random effects covariance matrix, and \eqn{\gamma} is a vector of two power parameters; see Zou (2006) for the adaptive lasso, and Hui et al. (2016) for regularized PQL selection in GLMMs using on adaptive lasso type penalties. 
 
If \code{cov.groups} is supplied, this means that some of the fixed effects coefficients should be treated and penalized collectively as a group. The most common cases where this is used is when you have factor or categorical variables with more than two levels, or when you have polynomial terms that should be dealt with together. For instance, suppose you have a model matrix consisting of six columns, where first three columns correspond to separate covariates (including the intercept) and the last three columns all correspond to dummy variables created for a factor variable with four levels , e.g. soil moisture with levels dry, moderately moist, very moist, wet. The coefficients from the last three columns should then be penalized together, and so we can set \code{cov.groups = c(1,2,3,4,4,4)}.  
 
In doing so, the adaptive lasso weights for the grouped coefficients are then constructed differently. Following on from the example above, we have the fixed effect weight for soil moisture defined as

\deqn{w^F = \|\tilde{\beta}\|^{-\gamma_1},}

where \eqn{\| \cdot \|} corresponds to the L2-norm and \eqn{\tilde{\beta}} are the fixed effect coefficients belonging in the group (three in this case). When entered into the \code{rpql} function, an adaptive group lasso (Wang and Leng, 2008) is applied to these set of coefficients, such that they are all encouraged to be shrunk to zero at the same time.

Of course, after construction the adaptive lasso weights can be manually altered before entering into the main \code{rpql} function e.g., if one wants certain fixed and/or random effects to not be penalized.
}

\value{
A list containing the following elements
  \item{fixef}{Fixed effect coefficient estimates from \code{lme4.fit}.}
  \item{ranef}{A list of random effect predicted coefficients from \code{lme4.fit}.}
  \item{ran.cov}{A list of random effects covariance matrices from \code{lme4.fit}.}
  \item{cov.groups}{The argument \code{cov.groups}. Defaults to \code{NULL}.}
  \item{pen.weights}{A list of adaptive lasso weights constructed from \code{lme4.fit}. Contains elements \code{pen.weight$fixed} and \code{pen.weights$random}, which are the weights for the fixed and random effects respectively. Please see details above as to their construction.}
  } 


\references{
\itemize{
  \item Hui, F.K.C., Mueller, S., and Welsh, A.H. (2016). Joint Selection in Mixed Models using Regularized PQL. Journal of the American Statistical Association: accepted for publication.
  \item Wang, H., and Leng, C. (2008). A note on adaptive group lasso. Computational Statistics & Data Analysis, 52, 5277-5286.
  \item Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American statistical association, 101, 1418-1429.
  }
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}

\section{Warnings}{
  \itemize{
  \item In order to construct sensible starting values and weights, this function should really only be used when \code{lme4.fit} is a fit of the saturated GLMM, i.e. all fixed and random effects included. 
  }
    } 

\seealso{
\code{\link{rpql}} for fitting and performing model selection in GLMMs using regularized PQL, which may use the values obtained from \code{build.start.fit} for starting values and adaptive lasso weights.
}

\examples{

##################
## Example 1: Bernoulli GLMM with grouped covariates. 
## Independent cluster model with 50 clusters and equal cluster sizes of 10
## Nine covariates where the last covariate (soil type) is a factor with four levels
n <- 50; p <- 8; m <- 10
set.seed(123)
X <- data.frame(matrix(rnorm(n*m*p),n*m,p), soil=sample(1:4,size=m*n,replace=TRUE))
X$soil <- factor(X$soil)
X <- model.matrix(~ ., data = X)
colnames(X) <- paste("X",1:ncol(X),sep="")

Z <- X[,1:5] ## Random effects model matrix taken as first five columns
true_betas <- c(-0.1,1,-1,1,-1,1,-1,0,0,0,0,0) 
true_D <- matrix(0,ncol(Z),ncol(Z))
true_D[1:3,1:3] <- matrix(c(9,4.8,0.6,4.8,4,1,0.6,1,1),
	3,3,byrow=TRUE) ## 3 important random effects 

simy <- gendat.glmm(id = list(cluster = rep(1:n,each=m)), X = X, beta = true_betas, 
	Z = list(cluster = Z), D = list(cluster = true_D), family = binomial())

  
\dontrun{
library(lme4)
dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
fit_satlme4 <- glmer(y ~ X - 1 + (Z - 1 | cluster), data = dat, 
	family = "binomial")
fit_sat <- build.start.fit(fit_satlme4, id = simy$id, gamma = 2, 
	cov.groups = c(1:9,10,10,10)) 

new.fit <- rpql(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, lambda = 0.01, 
	pen.type = "adl", pen.weights = fit_sat$pen.weights,
	cov.groups = fit_sat$cov.groups, start = fit_sat, family = binomial())  
	}

}

