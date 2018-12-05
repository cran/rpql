\name{rpql}

\alias{rpql}

\alias{rpql.default}

\alias{print.rpql}

\title{Joint effects selection in GLMMs using regularized PQL.}

\description{
\code{rpql} offers fast joint selection of fixed and random effects in Generalized Linear Mixed Model (GLMMs) via regularization. The penalized quasi-likelihood (PQL) is used as a loss function, and penalties are added on to perform fixed and random effects selection. This method of joint selection in GLMMs, referred to regularized PQL, is fast compared to information criterion and hypothesis testing (Hui et al., 2016). 

Please note \code{rpql} is the core workshops function that performed regularized PQL on a single set of tuning parameters. \code{rpqlseq} is a wrapper to permit a sequence of tuning parameter values. The latter is often what users want to use.
}

\usage{
rpql(y, ...)
  
\method{rpql}{default}(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, 
  pen.type = "lasso", start = NULL, cov.groups = NULL, pen.weights = NULL, 
  hybrid.est = FALSE, offset = NULL, intercept = TRUE, save.data = FALSE, 
  control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, 
  scad.a = 3.7, mcp.gamma = 2, seed = NULL), ...)

\method{print}{rpql}(x, ...)
  }

\arguments{
  \item{y}{A vector of responses}

  \item{X}{A model matrix corresponding to the fixed effects. It should have the same number of rows as the length of \code{y}. An intercept column must be included if a fixed intercept is desired.}
  
  \item{Z}{A list with each element being a model matrix for a set of random effects. Each element of \code{Z} is referenced by a vector of IDs given by the corresponding element in the list \code{id}. Each model matrix (element of \code{Z}) should have the same number of rows as the length of \code{y}.}
  
  \item{id}{A list with each element being a vector of IDs that reference the model matrix in the corresponding element in the list \code{Z}. Each vector of IDs \emph{must} be integers (but not factors).}

  \item{x}{An object for class "rpql".}

  \item{family}{The distribution for the responses in GLMM. The argument must be applied as a object of class "family". Currently supported arguments include: \code{gaussian()}, \code{poisson()}, \code{binomial()}, \code{Gamma()}, \code{nb2()} for negative binomial, \code{LOGNO()} for log-normal, and \code{ZIP()} for zero-inflated Poisson.}

  \item{trial.size}{The trial size if \code{family = binomial()}. Either takes a single non-zero value or a vector of non-zero values with length the same as the number of rows in \code{X}. The latter allows for differing trial sizes across responses. Defaults to 1.} 

  \item{lambda}{A vector of length one or two specifying the tuning parameters used in regularized PQL. If two elements are supplied, then first and second elements are for the fixed and random effects penalty respectively. If one element, then it is applied to both penalties.}

  \item{pen.type}{A vector of one or two strings, specifying the penalty used for variable selection. If two elements are supplied, then first and second strings are the fixed and random effects penalty respectively. If one element, the same type of penalty is used. Currently supported argument include: "\code{lasso}" for standard lasso (Tibshirani, 1996), "\code{scad}" for SCAD penalty with \eqn{a} controlled by \code{control$scad.a} (Fan and Li, 2001), "\code{adl}" for adaptive lasso (Zou, 06), "\code{mcp}" for MC+ penalty with \eqn{\gamma} controlled by controlled by \code{control$mcp.gamma} (Zhang, 2010). If the adaptive lasso is used, then \code{pen.weights} must also be supplied. Defaults to standard lasso penalty for both fixed and random effects.}

  \item{start}{A list of starting values. It must contain the following elements: \code{start$fixef} as starting values for the fixed effect coefficients, \code{start$ranef} which is a list containing matrices of starting values for the random effects coefficients. It may also contain \code{start$D} which is a list of matrices to act as starting values for random effects covariance matrices.} 
  
  \item{cov.groups}{A vector specifying if the columns of \code{X} (including the intercept) should be regarded and therefore penalized in groups. For example, if one or more of the fixed effect covariates are factors, then \code{lme4} will automatically create dummy variables in the model matrix and estimate coefficients for each level, using one level as the reference. \code{cov.groups} is then used to identify all the coefficients that corresponds to that factor, such that all of these coefficients are penalized collectively as a group. Defaults to NULL, in which case it is assumed all coefficients should be treated independently. Please see the details and examples for more details.}
  
  \item{pen.weights}{A list containing up to two elements for additional (adaptive lasso) weights to be included for penalization. This must be supplied if \code{pen.type} has one or both elements set to "adl", otherwise it is optional. A weights equal to zero implies no penalization is applied to the parameter. The two elements in the list are as follows: for fixed effects, \code{pen.type$fixed} should be a vector with length equal to the number of columns in \code{X}. For random effects, \code{pen.weights$ran} should be a list of the same length as the list \code{Z}, where each element in that list is a vector with length equal to the number of columns in the corresponding element of the list \code{Z} (recall that each element of \code{Z} is a model matrix). Defaults to NULL, in which case there are no weights involved in the penalization.}

  \item{hybrid.est}{Should a hybrid estimation approach be used? That is, once model selection is performed using regularized PQL, should the submodel be re-estimated using the \code{lme4} package, if possible? Defaults to FALSE.}

  \item{offset}{This can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. It should be numeric vector of length equal to \code{y}. Defaults to NULL.}

  \item{intercept}{Is one of the columns of \code{X} an intercept term? This is used to indicate the presence of a fixed intercept in the model, which subsequently will NOT be penalized. Defaults to TRUE.}

  \item{save.data}{Should \code{y, X}, and \code{Z}, be saved as part of the output? Defaults to FALSE. The data is not saved by default in order to save memory.}

  \item{control}{A list controlling the finer details of the rPQL algorithm. These include: 
	\itemize{
		\item{\code{tol}: }{Tolerance value for convergence in the regularized PQL to be declared, where convergence is measured as the difference between the estimated parameters in successive iterations. Defaults to a value of 1e-4.}
		
		\item{\code{maxit}: }{The maximum number of update iterations for regularized PQL. Defaults to 100.}
		
		\item{\code{trace}: }{Should the update estimates of the fixed effect coefficients and that random effect covariance matrices be printed at each iteration? Defaults to FALSE.}
		
		\item{\code{restarts}: }{The number of restarts to try in case the algorithm diverges, i.e. the fixed effect coefficients and /or random effects covariance matrices "blow up". Defaults to a value of 5. Divergence is mostly likely to occur when you have count responses with some extremely large counts in there, in which regularized PQL can throw a hissy fit.}
		
		\item{\code{scad.a, mcp.gamma}: }{Controls the \eqn{a} and \eqn{\gamma} parameters in the SCAD and MC+ penalty respectively. Defaults to \eqn{a = 3.7} (Fan and Li, 2001) and \eqn{\gamma = 2} (Zhang, 2010) respectively. Please note these parameters are only in use when \code{pen.type} involves these penalties.}

		\item{\code{seed}: }{A seed that can be used if results need to be replicated. Defaults to NULL, in which case a random seed is used.}
		}
	}

  \item{...}{Not used.}
}

\details{
\bold{Intro}

Generalized Linear Mixed Models (GLMMs) are an extension of Generalized Linear Models (GLM, see the \code{glm} function) to include one or more sets of random effects. For \eqn{i = 1,\ldots,n}, where \eqn{n} is the length of \code{y}, we have 

\deqn{g(\mu_{i}) = \bm{x}^T_i \bm{\beta} + \bm{z}^T_{i1} \bm{b}_{i1} + \bm{z}^T_{i2} \bm{b}_{i2} + \ldots,}

where \eqn{g(\cdot)} is the link function, \eqn{\mu_i} is the mean of the distribution for observation \eqn{i}, \eqn{\bm{x}_i} is row \eqn{i} of the fixed effects model matrix \code{X}, and \eqn{\bm{\beta}} is the fixed effects coefficients. For the random effects, \eqn{\bm{z}_{i1}} is row \eqn{i} of the random effects model matrix in the first element of \code{Z}, \eqn{\bm{z}_{i2}} is from the second element of \code{Z} and so forth. The random effects \eqn{\bm{b}_{i1}}, \eqn{\bm{b}_{i2}, \ldots} are drawn from a multivariate normal distribution with mean zero and differing covariance matrices \eqn{\bm{D}_1, \bm{D}_2, \ldots}.

Note that having lists for \code{id, Z}, allows for multiple sets of random effects to be included in the GLMM. This is analogous to the \code{lme4} package, where multiple random effects are permitted in the formula e.g., \code{(1|creek) + (1|creek:sample)}. If the GLMM contains only one set of random effects, e.g., in longitudinal data, then the two lists will all contain only one element. Cases where multiple sets of random effects may be used include nested and crossed designs, in which case \code{id, Z}, will have two or more elements. It is recommended that the user think through and design these lists carefully to ensure that they are actually constructing the appropriate GLMM of interest. Yes it takes some getting use too, and we apologize for this =( Please see examples below for some ideas.

%While they are a flexible class of models for handling non-normal, correlated data (see for instance Bolker et al., 2009, and Verbeke and Molenberghs, 2009), estimating and performing variable selection in GLMMs presents some major challenges. In the former, the marginal log-likelihood requires integration over the unobserved random effects, and this is computationally challenging. The \code{lme4} package, for instance, typically uses the Laplace's approximation to approximate the integral. In the latter, one may be interested in performing joint selection over the fixed and random effects selection. Methods like hypothesis testing using \code{anova} and information criterion may become computationally burdensome in such case, and also may not necessarily perform that well (see Mueller et al., 2013, for a review of variable selection in linear mixed models). 

\bold{Regularized PQL}

Regularized PQL is designed as a fast approach to joint selection to GLMMs (Hui et al., 2016). It works by taking the penalized quasi-likelihood (PQL, Breslow and Clayton, 1993) and adding on penalties to perform selection of the fixed and random effects. That is, maximize the regularized PQL function

\deqn{\ell = \sum\limits_{i=1}^n \log(f(y_i | \bm{\beta}, \bm{b}_{i1}, \bm{b}_{i2}, \ldots)) - \frac{1}{2} \sum\limits_{i=1}^n \bm{b}^T_{i1}\bm{D}^{-1}_1 \bm{b}_{i1} - \frac{1}{2} \sum\limits_{i=1}^n \bm{b}^T_{i2}\bm{D}^{-1}_2 \bm{b}_{i2} - \ldots - P_{\lambda}}

where \eqn{P_{\lambda}} denotes penalties to shrink the fixed effect \eqn{\bm{\beta}} and random effect \eqn{\bm{b}_{i1}}, \eqn{\bm{b}_{i2}, \ldots} coefficients, which depend on one or more tuning parameters \eqn{\lambda}. Like the PQL itself, regularized PQL is a fast approach for estimating GLMMs because it treats the random effects as "fixed" coefficients, and therefore no integration is required. Penalties are then used to shrunk one or more \eqn{\bm{\beta}}'s and \eqn{\bm{b}}'s to zero, the latter done so in a group-based manner, in order to perform joint selection (see Hui et al., 2016, for details). In short, regularized PQL is able to fit many GLMMs in a relatively short period of time, which in turn facilitates the construction of a solution or regularization path ranging from the null (intercept-only) to the full (saturated) model. A tuning parameter selection method such as information criterion can then be used to pick the select the final subset of fixed and random effects. A few penalty types are available in the package, from which we prefer to use the adaptive LASSO (with weights based on the full model, Zou, 2006) mainly because by having weights, we can avoids have to search through a two-dimensional grid of tuning parameter values.

Note that if one only wanted to penalize the fixed effects and leave the random effects unpenalized, this can be achieved by setting the second element/s of lambda equal to to e.g., \code{lambda = c(1,0)}. Note though that in longitudinal studies, for covariates included as both fixed and random effects, if the random effects is not penalized then neither should the fixed effect. This ensures that no covariates end up being selected in the model as a purely random effects (non-hierarchical shrinkage, Hui et al., 2016). This can be accounted for also setting the corresponding elements of \code{pen.weights$fixed} to zero.

\bold{AN IMPORTANT NOTE} 

While regularized PQL is relatively fast, it will produce biased estimates of the fixed and random effects parameters for non-normal responses, especially if the amount of data to estimate each random effect is not large e.g., if the number of time points or cluster size is not large. We envision regularized PQL as a method of joint variable selection ONLY, and strongly encourage the user to adopt a hybrid estimation approach (using \code{hybrid.est = TRUE}, for instance). That is, once model selection is performed using regularized PQL, the final submodel should be re-estimated using more exact methods like quadrature or MCMC.

Because regularized PQL treats the random effects as ``fixed" coefficients and therefore penalizes these, then the random effects covariance matrices \eqn{\bm{D}_1, \bm{D}_2, \ldots} are regarded more as nuisance parameters. This is in contrast to traditional maximum likelihood estimation where the random effect coefficients \eqn{\bm{b}_{i1}}, \eqn{\bm{b}_{i2}, \ldots} are integrated over. As nuisance parameters, regularized PQL employs an iterative estimator based on maximizing the Laplace-approximated marginal log-likelihood, assuming all other parameters are fixed, for estimating the covariance matrix \eqn{\bm{D}_1, \bm{D}_2, \ldots}. This iterative estimator was used in Hui et al., (2016) for independent clustered data specifically. When they are multiple sets of random effects, each covariance matrix is estimated conditionally on all others i.e., the random effect coefficients corresponding to all other random effects are held constant. This can be thought of as employing a series of conditional Laplace approximations to obtain updates for \eqn{\bm{D}_1, \bm{D}_2, \ldots}.

%As nuisance parameters, regularized PQL currently offers three possible estimates of the \eqn{\bm{D}_1, \bm{D}_2, \ldots}. The default is a simple sample covariance estimator based on the estimated random effects coefficients, "\code{bb}". That is, given estimates of \eqn{\bm{b}_{i1}} for all \eqn{i = 1,\dots,n} based on regularized PQL, the estimate of \eqn{\bm{D}_1} is given by 

%\deqn{\hat{\bm{D}}_1 = \frac{1}{n} \sum\limits_{i=1}^n \hat{\bm{b}}_{i1} \hat{\bm{b}}^T_{i1},} 

%and similarly for \eqn{\bm{D}_2, \ldots}. The estimate of the random effects covariance matrix is consistent, provided the number of observations available to estimate each of the \eqn{\bm{b}_{i1}}'s i.e., the cluster size, grows with \eqn{n}. The two other options for estimating the covariance matrix are: 1) a weighted sample covariance estimator, where the weights are based on the number of observations available to estimate each random effect coefficient i.e., cluster size, such that more weight is given to coefficients estimated with a larger cluster size. This estimator is useful when you have unbalanced cluster sizes, which is almost always the case in real life, and 2) an iterative estimator based on maximizing the Laplace-approximated marginal log-likelihood assuming all other parameters are fixed. This iterative estimator was used in Hui et al., (2016), but is only applicable when you have independent clustered data e.g., longitudinal studies.

\bold{A not so short discussion about information criterion}

How to choose the tuning parameters for penalized regression is an active area of area of research in statistics (see for instance Zhang et al., 2010, Hui et al., 2014), with the most popular solutions being cross validation and information criteria. That is, a solution path is constructed and the best submodel is then chosen by minimizing the value of the information criterion. Anyway, \code{rpql} offers the following information criteria for tuning parameter selection, as available in \code{ics} in the output. Please note all of the criteria below use only the first part of the PQL function as the loss function i.e., \eqn{IC = -2\sum\limits_{i=1}^n \log(f(y_i | \bm{\beta}, \bm{b}_{i1}, \bm{b}_{i2}, \ldots)) +} model complexity terms. 

\enumerate{
\item A AIC-type criterion that penalizes a values of 2 for every non-zero fixed effect coefficient, and, for each set of random effects, penalizes a value of 2 for every non-zero random effect coefficient in that set.

\item A BIC-type criterion that penalizes a value of \eqn{\log(n)} for every non-zero fixed effect coefficient, and, for each set of random effects, penalizes a value of \eqn{\log(n_c)} for every non-zero, unique element in covariance matrix for that set, where \code{n_c} denotes the number of clusters corresponding to that random effect.

\item A BIC-type criterion that penalizes a value of \eqn{\log(n)} for every non-zero fixed effect coefficient, and, for each set of random effects, penalizes a value of \eqn{\log(n)} for every non-zero, unique element in covariance matrix for that set. This combination of penalties is the one used in the package \code{lme4}. 

\item Three hybrid information criteria that penalizes a value \eqn{\log(n)} for every non-zero fixed effect coefficient, and, for each set of random effects, penalizes a value of 2/1/0.5 for every non-zero random effect coefficient in that set.
}

Selection consistency for all but the first AIC criteria have been established, although empirically performance may differ. We generally prefer the three hybrid criterion, although it is recommended that the user tries several of them and see how results differ! 

}

\value{
An object of class "rpql" containing the following elements:
\item{call}{The matched call.}

\item{fixef}{A vector of estimated fixed effect coefficients, \eqn{\beta}.}

\item{ranef}{A list with each element being a matrix of estimated (predicted) random effect coefficients, \eqn{\bm{b}_{i1}}, \eqn{\bm{b}_{i2}}, and so on.}

\item{ran.cov}{A list with each element being an estimated random effect covariance matrices, \eqn{\bm{D}_1, \bm{D}_2, \ldots}.}

\item{logLik}{The (unpenalized) PQL likelihood value at convergence.}

\item{phi, shape, zeroprob}{Estimates of nuisance parameters (if appropriate), including the variance and overdispersion parameter for normal, lognormal and negative binomial families, the shape parameter for the Gamma family, and the probability of a structural zero for zero-inflated Poisson family.}

\item{family}{The family fitted.}

\item{n}{The length of \code{y}.}

\item{id}{The \code{id} argument.}

\item{lambda, pen.type}{The tuning parameters and penalties used.}

\item{ics}{A vector containing the number of estimated parameters in the GLMM (note regularized PQL treats the random effects as "fixed"), and some information criteria. Please see \code{details} above for more information.}

\item{nonzero.fixef}{A vector indexing which of the estimated fixed effect coefficients are non-zero.}

\item{nonzero.ranef}{A list with each element being a vector indexing which of the estimated random effects are non-zero, i.e. which of the diagonal elements in the corresponding element of \code{ran.cov} are non-zero.}

\item{hybrid}{The estimated fit from \code{lme4}, if \code{hybrid.est = TRUE}.}

\item{y,X,Z}{The data the GLMM is fitted to, if \code{save.data = TRUE}.}
}

\references{
\itemize{
   \item Breslow, N. E., and Clayton, D. G. (1993). Approximate inference in generalized linear mixed models. Journal of the American Statistical Association, 88, 9-25.
   
   \item Fan, J., and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. Journal of the American statistical Association, 96, 1348-1360.
   
   \item Hui, F.K.C., Mueller, S., and Welsh, A.H. (2017). Joint Selection in Mixed Models using Regularized PQL. Journal of the American Statistical Association, 112, 1323-1333.
   
   \item Hui, F.K.C., Mueller, S., and Welsh, A.H. (2017). Hierarchical Selection of Fixed and Random Effects in Generalized Linear Mixed Models. Statistica Sinica, 27, 501-518.

   \item Hui, F. K., Warton, D. I., and Foster, S. D. (2014). Tuning parameter selection for the adaptive lasso using ERIC. Journal of the American Statistical Association, 110, 262-269.   
	
   \item Lin, X., and Breslow, N. E. (1996). Bias correction in generalized linear mixed models with multiple components of dispersion. Journal of the American Statistical Association, 91, 1007-1016.
	
   \item Mueller, S., Scealy, J. L., and Welsh, A. H. (2013). Model selection in linear mixed models. Statistical Science, 28, 135-167.
   
   \item Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society. Series B (Methodological), 58, 267-288.
   
   \item Zhang, Y., Li, R., and Tsai, C. L. (2010). Regularization parameter selections via generalized information criterion. Journal of the American Statistical Association, 105, 312-323.
   
   \item Zhang, C. H. (2010). Nearly unbiased variable selection under minimax concave penalty. The Annals of Statistics, 38, 894-942.
   
   \item Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American statistical association, 101, 1418-1429.
  }
}

\author{
\packageAuthor{rpql}

Maintainer: \packageMaintainer{rpql}
}

\section{Warnings}{
\itemize{
\item We strongly recommend you scale your responses (if normally distributed) and any continuous covariates, otherwise \code{rpql} like all penalized likelihood methods, may not make much sense!

\item Like its standard unpenalized counterpart, regularized PQL can produce very bias parameter estimates in finite samples, especially if you do not have a lot of data to estimate each random effect. We therefore envision regularized PQL as a tool for fast model selection in GLMMs, and strongly recommend you re-estimate the final submodel using more accurate estimation methods i.e., use a hybrid estimation approach, in order to obtain better final parameter estimates and predictions of the random effects.

\item If \code{save.data = TRUE}, the data you fitted the GLMM is also saved as part of the output, and this can potentially take up a lot of memory. 

\item If you are constantly suffering convergence issues with regularized PQL, even after multiple restarts, consider increasing \code{lambda[2]} to penalized the random effects more and stabilize the estimation algorithm. You may also want to consider better starting values, in particular, smaller values of \code{start$ranef}. Good luck!

}
}

\seealso{
\code{\link{rpqlseq}} for the wrapper function that runs \code{rpql} multiple times on a sequence of tuning parameter values, \code{\link{build.start.fit}} for building \code{start} lists from a GLMM fitted using the \code{lme4} package, \code{\link{summary}} for a summary of the regularized PQL fit. For alternative methods of fitting GLMMs, you may also want be check out the packages \code{lme4, nlme, MCMCglmm} and \code{glmmADMB}.
}

\examples{
## Please note all examples below use the \code{rpqlseq} wrapper function. 

library(lme4)
library(gamlss.dist)

##################
## Example 1: Poisson GLMM on simulated data 
## Indepenent cluster model with 30 clusters and equal cluster sizes of 10
## 9 fixed and random effect covariates including a fixed and random intercept
library(mvtnorm)
set.seed(1)
n <- 30; m <- 10; p <- 8; 
## Generate rows of a model matrix from a multivariate normal distribution 
## with AR1 covariance structure. 

H <- abs(outer(1:p, 1:p, "-")) 
X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 
Z <- X 
true_betas <- c(0.1,1,-1,-1,1,rep(0,p-4)) ## 5 truly important fixed effects
true_D <- matrix(0,ncol(Z),ncol(Z))
true_D[1:3,1:3] <- matrix(c(1,0.6,0.6,0.6,1,0.4,0.6,0.4,1),3,3,byrow=TRUE) 
## 3 important random effects

simy <- gendat.glmm(id = list(cluster=rep(1:n,each=m)), X = X, beta = true_betas, 
	Z = list(cluster=Z), D = list(cluster=true_D), family = poisson()) 

	
\dontrun{
## Construct a solution path using adaptive LASSO for selection 
dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
fit_satlme4 <- glmer(y ~ X - 1 + (Z - 1 | cluster), data = dat,
	family = "poisson")
fit_sat <- build.start.fit(fit_satlme4, gamma = 2)
## Please see example 3 for another way of constructing the adaptive weights

lambda_seq <- lseq(1e-6,1,length=100)
fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
	family = poisson(), lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)

summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs

## Note, if you wanted to penalized the fixed effects only, this can achieved
## by setting fit_sat$pen.weights$random$cluster <- rep(0,ncol(simy$Z$cluster))

## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
XMM <- unname(model.matrix(fit_satlme4)) 
ZMM <- getME(fit_satlme4,"mmList"); names(ZMM) <- "cluster"
lambda_seq <- lseq(1e-6,1,length=100)
fit <- rpqlseq(y = simy$y, X = XMM, Z = ZMM, id = simy$id, 
 	family = poisson(), lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)
## Big thanks for Andrew Olney for this suggestion!
}


##################
## Example 2: Similar to example 1 but with Bernoulli GLMMs 
## 30 clusters, cluster size of 20
library(mvtnorm)
set.seed(1)
n <- 30; m <- 20; p <- 8; 
## Generate rows of a model matrix from a multivariate normal distribution 
## with AR1 covariance structure. 

H <- abs(outer(1:p, 1:p, "-")) 
X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 
Z <- X 
true_betas <- c(-0.1,1,-1,1,-1,rep(0,p-4)) ## 5 truly important fixed effects
true_D <- matrix(0,ncol(Z),ncol(Z))
true_D[1:3,1:3] <- diag(c(3,2,1), nrow = 3)
## 3 important random effects

simy <- gendat.glmm(id = list(cluster=rep(1:n,each=m)), X = X, 
  beta = true_betas, Z = list(cluster=Z), D = list(cluster=true_D), family = binomial()) 

	
\dontrun{
## Construct a solution path using adaptive LASSO for selection 
dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
fit_satlme4 <- glmer(y ~ X - 1 + (Z - 1 | cluster), data = dat, 
	family = "binomial")
fit_sat <- build.start.fit(fit_satlme4, gamma = 2)

lambda_seq <- lseq(1e-6,1,length=100)
best.fit <- list(ics = rep(Inf,6))
fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)
	
summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs

## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
XMM <- unname(model.matrix(fit_satlme4)) 
ZMM <- getME(fit_satlme4,"mmList"); names(ZMM) <- "cluster"
lambda_seq <- lseq(1e-6,1,length=100)
fit <- rpqlseq(y = simy$y, X = XMM, Z = ZMM, id = simy$id, 
 	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)
}


##################
## Example 3: Bernoulli GLMMs on simulated data
## Nested data with 200 observations in total: split into 10 creeks, 
## 5 samples nested within each creek
## Please see example in gendat.glmm for further details
mn <- 200; 
X <- matrix(1,mn,1); 
ids <- list(samples = rep(1:50,each=4), creek = rep(1:10,each=20)) 
## We have two sets of random intercepts only, one for creek and one for 
## samples nested within creek.
Zs <- list(samples = X, creek = X) 

true_betas <- -0.1 
true_D <- list(samples = as.matrix(0.001), creek = as.matrix(1)) 
## Please ensure each element of true_D is a matrix

simy <- gendat.glmm(id = ids, X = X, beta = true_betas, Z = Zs, 
	D = true_D, trial.size = 1, family = binomial())

\dontrun{
## Construct a solution path use adaptive LASSO for selection
## Here is another way of constructing the adaptive weights:
## Use the fact that rpql can do a final fit based on maximum likelihood
## to obtain a good saturated fit.
fit_sat <- rpql(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
	family = binomial(), lambda = 0, hybrid = TRUE)
fit_sat <- build.start.fit(fit_sat$hybrid, gamma = 2)
	
fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)

summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
}
   

##################
## Example 4: Linear mixed models on Alfalfa split-plot data

\dontrun{

library(nlme)
data(Alfalfa)
Alfalfa$Yield <- scale(Alfalfa$Yield)
X <- as.matrix(model.matrix(~ Date, data = Alfalfa)) 
## Note Date is categorical variable!
colnames(X)[1] <- "x1"
Z <- list(BlockVariety = matrix(1,nrow(X),1), Block = matrix(1,nrow(X),1))
## Four samples of each Block*Variety
ids <- list(BlockVariety = rep(1:(nrow(X)/4),each=4), 
	Block = as.numeric(Alfalfa$Block)) 

## How you would fit it in lme4
fit_satlme4 <- lmer(Yield ~ X - 1 + (1|Block/Variety), data = Alfalfa)
fit_sat <- build.start.fit(fit_satlme4, cov.groups = c(1,2,2,2), gamma = 2)

## Construct a solution path using adaptive LASSO for selection
lambda_seq <- lseq(1e-5,2,length=100)
fit <- rpqlseq(y = Alfalfa$Yield, X = X, Z = Z, id = ids, 
	lambda = lambda_seq, cov.groups = c(1,2,2,2), pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)

summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs


## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
X <- unname(model.matrix(fit_satlme4)) 
Z <- getME(fit_satlme4, "mmList"); names(Z) <- c("BlockVariety", "Block")
lambda_seq <- lseq(1e-6,1,length=100)
fit <- rpqlseq(y = Alfalfa$Yield, X = X, Z = Z, id = ids, 
	lambda = lambda_seq, cov.groups = c(1,2,2,2), pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)
}


##################
## Example 5: Linear mixed models on sleep study dataset

\dontrun{
data(sleepstudy)

## How you fit it in lme4
## Response is scaled so as to avoid large variances and easier intepretation
sleepstudy$Reaction <- scale(sleepstudy$Reaction) 
sleepstudy$Days <- scale(sleepstudy$Days)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

## How you fit it using rpql
## Construct a solution path using adaptive LASSO for selection 
X <- cbind(1,sleepstudy$Days)
Z <- list(subject = X)
ids <- list(subject = as.numeric(sleepstudy$Subject))
fit_sat <- build.start.fit(fm1, gamma = 2)

lambda_seq <- lseq(1e-4,1,length=100)
fit <- rpqlseq(y = sleepstudy$Reaction, X = X, Z = Z, id = ids, 
	lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)

summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
## Best fit might well be the saturated fit! 
## This is at least consistent with confint(fm1)

## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
X <- unname(model.matrix(fm1)) 
Z <- getME(fm1, "mmList"); names(Z) <- "subject"
lambda_seq <- lseq(1e-6,1,length=100)
fit <- rpqlseq(y = sleepstudy$Reaction, X = X, Z = Z, id = ids, 
	lambda = lambda_seq, pen.type = "adl", 
	pen.weights = fit_sat$pen.weights, start = fit_sat)
}


##################
## Example 6: GLMM with lognormal responses
## Fixed effects selection only

\dontrun{
n <- 50; m <- 10; p <- 8; 
H <- abs(outer(1:p, 1:p, "-")) 
X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 
Z <- X[,1:3] ## 3 random effects all of which important
true_betas <- c(0.1,1,-1,-1,1,rep(0,p-4)) ## 5 important fixed effects
true_D <- matrix(0,ncol(Z),ncol(Z))
true_D[1:3,1:3] <- matrix(c(1,0.6,0.6,0.6,1,0.4,0.6,0.4,1),3,3,byrow=TRUE) 

simy <- gendat.glmm(id = list(cluster=rep(1:n,each=m)), X = X, 
	beta = true_betas, Z = list(cluster=Z), D = list(cluster=true_D), 
	family = LOGNO(), phi = 1) 

## We will use the lasso penalty for fixed effects only with no weights
## Note lognormal mixed models are usually hard to fit by maximum likelihood in R!
## Hence adaptive weights are sightly hard to obtain

## Note also that since random effects are not penalized, then generally 
## the corresponding fixed effect covariates should not be penalized 
## (at least in longitudinal studies), in keeping in line with the 
## hierarchical principle of the effects.
## To account for this in the above, we can use the pen.weights argument 
## to prevent penalization of the first three fixed effect covariates

fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
  family = LOGNO(), lambda = lambda_seq, pen.type = "lasso", start = NULL, 
	pen.weights = list(fixed = rep(c(0,1), c(3,ncol(X)-3))))

summary(fit$best.fit[[3]])  
# apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
}

}
