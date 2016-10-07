############
## Regularized PQL (with a penalty on fixed effects and a group penalty on the random effects
## Application to independent cluster models; i = 1,..,n; j = 1,...m
## g(\eta_i) = \bm{x}'_i \bm{\beta} + \bm{z}'_i \bm{b}_i; 
## \sum\limits_{i=1,j=1}^{n,m} log f(y_{ij}|\bm{b}_i, \bm{\beta}) - 0.5\sum\limits_{i=1}^n \bm{b}'_i D \bm{b}_i - \lambda \sum\limits_{k=1}^{p_f} |\beta_k| - \lambda \sum\limits_{k=1}^{p_r} ||\bm{b}'_k \bm{b}_k ||^{1/2}

## Version 5
## optimization details nested in a control command, and max.iter renamed maxit
## Default # of restarts lowered to 5
## Scad a parameter and MC+ gamma parameter now can be customized within the control argument
## rpqlseq now created, which wraps the rpql default function to permit a sequence of tuning parameter values

## TODO: 1) consider using threshold operators, including LLA and threshold; 2) rpqlseq to run a hybrid.est at the end if possible; 2) create a simulate function that runs gendat.glmm on a fitted rpql model?
##############
# library(mvtnorm)
# library(Matrix)
# library(MASS)
# library(lme4)
# library(gamlss.dist)
# source("auxilaryfunctions.R")



## pen.weights is a list with up to two elements for adl; first contains weights for fixed, and second contains a list of weights for random. 
# y = sim.y$y; X = sim.y$X; Z = sim.y$Z; id = sim.y$id; family = binomial(); lambda = 0.1; pen.type = "mcp"; hybrid.est = FALSE; offset = NULL; trial.size = 1; pen.weights = NULL; cov.groups = NULL; trace = TRUE; control$restarts <- 10; control$seed = NULL; intercept = TRUE; save.data = FALSE; tol = 1e-4; control$maxit = 100; start = NULL

rpql <- function(y, ...) UseMethod("rpql")


rpqlseq <- function(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, pen.type = "lasso", start = NULL, cov.groups = NULL, pen.weights = NULL, offset = NULL, intercept = TRUE, save.data = FALSE, control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, seed = NULL), ...) {
	lambda <- as.matrix(lambda)
	if(ncol(lambda) == 1) lambda[,1] <- sort(lambda[,1], decreasing = FALSE) 
	
	if(!("tol" %in% names(control))) control$tol <- 1e-4
	if(!("maxit" %in% names(control))) control$maxit <- 100
	if(!("trace" %in% names(control))) control$trace <- FALSE
	if(!("restarts" %in% names(control))) control$restarts <- 5
	if(!("scad.a" %in% names(control))) control$scad.a <- 3.7
	if(!("mcp.gamma" %in% names(control))) control$mcp.gamma <- 2
	if(!("seed" %in% names(control))) control$seed <- NULL 
	
	## Do a starting fit for formatting reasons
	start.fit <- rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[1,], pen.type = "lasso", start = start) 

	collect.ics <- matrix(0, nrow = nrow(lambda), ncol = length(start.fit$ics)+1) ## First column in start.fit$ics is DoF
	colnames(collect.ics) <- c("PQL Likelihood", names(start.fit$ics)) 
	best.fits <- vector("list", length(start.fit$ics)-1); 
	for(k in 1:length(best.fits)) best.fits[[k]] <- start.fit
	names(best.fits) <- names(start.fit$ics)[-1]
	rm(start.fit)		
			
	for(l1 in 1:nrow(lambda)) {
		message("Onto ", l1)
		if(l1 == 1) prev.fit <- start
		new.fit <- rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[l1,], pen.type = pen.type, pen.weights = pen.weights, start = prev.fit, cov.groups = cov.groups, hybrid.est = FALSE, offset = offset, intercept = intercept, save.data = save.data, control = control)
		
		for(l2 in 1:length(best.fits)) { if(new.fit$ics[l2+1] < best.fits[[l2]]$ics[l2+1]) best.fits[[l2]] <- new.fit }
		prev.fit <- new.fit
		collect.ics[l1,] <- c(new.fit$logLik1, new.fit$ics)
		}
	
	out <- list(best.fits = best.fits, collect.ics = collect.ics, lambda = lambda)	
	}


rpql.default <- function(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, pen.type = "lasso", start = NULL, cov.groups = NULL, pen.weights = NULL, hybrid.est = FALSE, offset = NULL, intercept = TRUE, save.data = FALSE, control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, seed = NULL), ...) {
	conv.eps <- 1e-3; #ran.cov.est <- "laplace"
	if(!("tol" %in% names(control))) control$tol <- 1e-4
	if(!("maxit" %in% names(control))) control$maxit <- 100
	if(!("trace" %in% names(control))) control$trace <- FALSE
	if(!("restarts" %in% names(control))) control$restarts <- 5
	if(!("scad.a" %in% names(control))) control$scad.a <- 3.7
	if(!("mcp.gamma" %in% names(control))) control$mcp.gamma <- 2
	if(!("seed" %in% names(control))) control$seed <- NULL 
	
	
	y <- as.vector(y); X <- as.matrix(X); num.fixed <- ncol(X); 
	if(det(t(X)%*%X) == 0) 
		warnings("Is X rank-deficient? Please double check!") 
	if(!is.list(Z) || !is.list(id)) 
		stop("Please supply id and Z as lists. The element in the two lists correspond to a vector of IDs and the random effects model matrix for those IDs respectively.")
	if(sum(duplicated(names(id))) > 0 || sum(duplicated(names(Z))) > 0) 
		stop("Please ensure the names of the elements in the lists id (and Z) are unique.")
	if(length(Z) != length(id)) 
		stop("The number of elements in Z and id should be the same. Thanks.")
# 	if(any(apply(X,2,sd) != 1) || any(apply(X,2,sd) != 1))
# 		message("You may want to consider standardizing your X and Z matrix prior to penalization.")
		
		
	if(is.null(offset)) offset <- rep(0,length(y))
	if(length(offset) != length(y)) stop("Offset should have the same length as y. Thanks.")

	
	num.ran <- get.lengths.id <- get.nrows.Z <- length.uniids <- numeric(length(Z)); 
	for(k in 1:length(Z)) { 
		Z[[k]] <- as.matrix(Z[[k]]); num.ran[k] <- ncol(Z[[k]]); get.nrows.Z[k] <- nrow(Z[[k]]) 
		}
	for(k in 1:length(id)) { 
		id[[k]] <- as.integer(id[[k]]); get.lengths.id[k] <- length(id[[k]]); length.uniids[k] <- length(unique(id[[k]])) 
		}
	if((get.nrows.Z != get.lengths.id) || (get.nrows.Z != nrow(X)) || (nrow(X) != get.lengths.id)) 
		stop("The number of rows in X, the length of each element in list id, and the number of rows in each element in list Z should all be the same. Thanks.")
	rm(get.lengths.id, get.nrows.Z)

	
	if(!is.null(cov.groups)) {
		actual.cov.groups <- cov.groups 
		if(control$trace) 
			message("Group penalization performed as cov.groups has been supplied.")
		if((ncol(X) != length(cov.groups))) 
			stop("If supplied, the length of cov.groups should be equal to the number of columns in X. Thanks") 
		}
	if(is.null(cov.groups)) { actual.cov.groups <- 1:ncol(X) }
		
		
	if(length(pen.type) == 1) pen.type <- rep(pen.type,2)
	if(sum(pen.type %in% c("lasso","scad","adl","mcp")) != 2) 
		stop("Current version does not permit specified penalty. Sorry!")		
	if(pen.type[1] == "adl" & !is.vector(pen.weights$fixed)) 
		stop("Please supply weights for the fixed effects in adaptive lasso. Specifically, pen.weights$fixed should be a vector.")
	if(pen.type[2] == "adl" & !is.list(pen.weights$ran)) 
		stop("Please supply weights for the random effects in adaptive lasso. Specifically, pen.weights$ran should be a list with the same number of elements as lists Z and id.")
	if(!is.null(pen.weights$fixed) & length(pen.weights$fixed) != ncol(X)) 
		stop("Weights provided for fixed effects must have the same length as the number of columns in X. Thanks.")
	if(!is.null(pen.weights$ran)) { 
		get.lengths.weights <- sapply(pen.weights$ran, length)
		if(!all(get.lengths.weights == num.ran)) 
			stop("The length of each element in pen.weights$ran should equal the number of columns in the corresponding element of list Z. Thanks.") 
		rm(get.lengths.weights) 
		}
	
	
# 	if(!(ran.cov.est %in% c("bb","wbb","laplace"))) 
# 		stop("Current version does not permit specified estimator of random effects covariance matrix. Sorry!")
# 	if(ran.cov.est == "laplace" & length(Z) > 1) {
# 		message("Laplace type estimator of random effects covariance matrix is permitted only if there is one set of random effects, i.e. length(Z) == length(id) = 1.") 
# 		message("Switching to ran.cov.est = wbb") 
# 		}
	
	
	if(!(family$family[1] %in% c("gaussian","poisson","binomial","negative.binomial","Gamma","LOGNO","ZIP"))) 
		stop("Current version does not permit specified family. Sorry!")
	if(family$family[1] == "binomial" & !(length(trial.size) %in% c(1,nrow(X)))) 
		stop("Binomial family requires trial.size to either be a single non-zero value or a vector of non-zero values with length equal to the number of rows in X.")
	if(family$family[1] != "binomial" & length(trial.size) != 1) {
		trial.size <- 1
		message("trial.size argument ignored for non-binomial families.")
		}
	
	
	if(length(lambda) == 1) lambda <- rep(lambda,2); lambda <- as.vector(unlist(lambda))
	if(!(length(lambda) %in% c(1,2))) stop("lambda should be a vector of length 1 or 2.")
	old <- .Random.seed
	on.exit( { .Random.seed <<- old } )	
	if(!is.null(control$seed)) set.seed(control$seed)
		
	if(is.null(colnames(X))) colnames(X) <- paste("x",1:ncol(X),sep=""); 
	for(k in 1:length(Z)) { colnames(Z[[k]]) <- paste("z",k,1:ncol(Z[[k]]),sep="") }
	
	
	init.fit.run <- 0
	if(is.null(start)) {
		get.init <- start.fixed(y = y, X = X, Z = Z, id = id, family = family, offset = offset, trial.size = trial.size)
		start <- build.start.fit.rpql(fit = get.init, id = id, num.ran = num.ran, cov.groups = cov.groups)
		init.fit.run <- 1
		}

		
	new.beta <- beta <- start$fixef
	new.b <- b <- start$ranef
	if(!is.list(new.b)) 
		stop("start$ranef should be a list of random effects. Thanks.")
	for(k in 1:length(new.b)) { new.b[[k]] <- b[[k]] <- as.matrix(b[[k]]) }
	new.phi <- phi <- new.shape <- shape <- new.zeroprob <- zeroprob <- 1
	if(family$family[1] == "gaussian") new.phi <- phi <- mean((y-X%*%beta-offset)^2) 
	if(family$family[1] == "LOGNO") new.phi <- phi <- mean((log(y)-X%*%beta-offset)^2) 
	if(family$family[1] == "negative.binomial") new.phi <- phi <- 0.01 
	if(family$family[1] == "Gamma") new.shape <- shape <- 1
	if(family$family[1] == "ZIP") new.zeroprob <- zeroprob <- 0.1

	if(is.null(start$ran.cov)) { 
		new.D <- D <- vector("list",length(id))
		for(k in 1:length(D)) new.D[[k]] <- D[[k]] <- cov(new.b[[k]])+diag(x=rep(1e-4,num.ran[k])) 
		}
	if(!is.null(start$ran.cov)) { new.D <- D <- start$ran.cov }
	

	## Construcs big block diagonal matrix based on an a vector of id and a Z matrix. Basically converts Z to a block Z.
	make.bigZ <- function(Z, id) {
		matlist <- split(as.data.frame(Z), id)
		matlist <- lapply(matlist, function(x) as.matrix(x))
		bigZ <- matrix(0, nrow = nrow(Z), ncol = ncol(Z)*length(matlist))
		for(k in 1:length(matlist)) { bigZ[as.numeric(rownames(matlist[[k]])), (k*ncol(Z)-ncol(Z)+1):(k*ncol(Z))] <- matlist[[k]] }
		bigZ <- Matrix(bigZ, sparse = TRUE)
		return(bigZ)
		}

		
	diff <- 1e4; cw.rpql.loglik <- 0; counter <- true.counter <- restart.counter <- 0
	ranef.fullshrunk <- rep(0,length(Z)) ## Indicator if random effects has been shrunk to zero. If 1, then ignore updating it!
	if(family$family[1] %in% c("LOGNO","gaussian")) { XtWX <- t(X)%*%X; } ## Define these outside as they do not change with iteration

	
	while(diff > control$tol & true.counter < control$maxit & restart.counter <= control$restarts) {
		all.Zb <- numeric(length(y)) 
		for(k2 in 1:length(id)) { all.Zb <- all.Zb + as.numeric(make.bigZ(Z = Z[[k2]], id = id[[k2]])%*%c(t(b[[k2]]))) }

		
		if(!(family$family[1] %in% c("gaussian","LOGNO"))) {
			new.etas <- X%*%beta + all.Zb + offset
			new.etas[new.etas > 30] <- 30; new.etas[new.etas < -30] <- -30

			if(family$family[1] == "negative.binomial") 
				new.weights <- (family$mu.eta(new.etas))^2/(family$variance(family$linkinv(new.etas),phi=phi)) ## 1/(variance*g'(mu)^2)
			if(family$family[1] == "Gamma")
				new.weights <- (family$mu.eta(new.etas))^2/(trial.size*family$variance(family$linkinv(new.etas))/shape) 
			if(family$family[1] == "ZIP") {
				new.weights <- (poisson()$mu.eta(new.etas))^2/(poisson()$variance(poisson()$linkinv(new.etas))) 
				weightmulti <- (1-zeroprob)*dpois(y,lambda=family$mu.linkinv(new.etas))/dZIP(x=y,mu=family$mu.linkinv(new.etas),sigma=zeroprob)
				weightmulti[!is.finite(weightmulti)] <- 0
				new.weights <- new.weights*weightmulti
				rm(weightmulti)
				}
			if(!(family$family[1] %in% c("ZIP","negative.binomial","Gamma"))) 
				new.weights <- (family$mu.eta(new.etas))^2/(trial.size*family$variance(family$linkinv(new.etas))) ## 1/(variance*g'(mu)^2)
				
			if(family$family[1] != "ZIP") 
				new.workresp <- new.etas + (y-trial.size*family$linkinv(new.etas))/family$mu.eta(new.etas) ## z = eta + (y-mu)*g'(mu)
			if(family$family[1] == "ZIP") 
				new.workresp <- new.etas + (y-family$mu.linkinv(new.etas))/poisson()$mu.eta(new.etas) 

			new.weights <- as.vector(new.weights)	
			new.weights[new.weights > 1e4] <- 1e4
# 			new.weights[new.weights < 1e-4] <- 1e-4
			XtWX <- t(X)%*%(X*as.vector(new.weights))
			XtWZ <- t(X)%*%((new.workresp-all.Zb-offset)*as.vector(new.weights))
			}

			
		## Update fixed effects conditional on random
		absbetas <- sapply(split(as.vector(beta),actual.cov.groups), l2.norm)[actual.cov.groups] + 1e-6
		D1 <- matrix(0,length(absbetas),length(absbetas))		
		if(pen.type[1] %in% c("adl","lasso")) diag(D1) <- length(y)*lambda[1]/as.vector(absbetas)
		if(pen.type[1] == "scad") diag(D1) <- length(y)*scad.deriv(absbetas, lambda=lambda[1], a = control$scad.a)/as.vector(absbetas)
		if(pen.type[1] == "mcp") diag(D1) <- length(y)*mcp.deriv(absbetas, lambda=lambda[1], gamma = control$mcp.gamma)/as.vector(absbetas)
		if(!is.null(pen.weights$fixed)) diag(D1) <- diag(D1)*c(pen.weights$fixed)
		if(intercept) diag(D1)[1] <- 0 
		D1[which(D1 > 1e8)] <- 1e8

		if(family$family[1] == "gaussian") 
			new.beta <- solve(XtWX + phi*D1)%*%t(X)%*%c(y-all.Zb-offset)
		if(family$family[1] == "LOGNO") 
			new.beta <- solve(XtWX + phi*D1)%*%t(X)%*%c(log(y)-all.Zb-offset)
		if(!(family$family[1] %in% c("gaussian","LOGNO"))) {
			new.beta <- try(solve(XtWX + D1)%*%XtWZ, silent=TRUE)
			if(inherits(new.beta,"try-error")) new.beta <- ginv(XtWX + D1)%*%XtWZ
			}
		if(lambda[1] > 0) {
			newbeta.l2 <- sapply(split(as.vector(new.beta),actual.cov.groups), l2.norm)
			get.zeros <- which(newbeta.l2 < conv.eps)
			new.beta[which(actual.cov.groups %in% get.zeros)] <- 0
			rm(get.zeros, newbeta.l2)
			}
			

		## Update random effects conditional on fixed
		for(k in 1:length(id)) {
			if(ranef.fullshrunk[k] == 1) next;
			
			all.Zb <- numeric(length(y)); 
			if(length(Z) > 1) { 
				for(k2 in (1:length(id))[-k]) all.Zb <- all.Zb + as.numeric(make.bigZ(Z = Z[[k2]], id = id[[k2]])%*%c(t(b[[k2]]))) 
				} ## Conditioning on all linear predictors for other ids, if is more than one id
			bigZ <- make.bigZ(Z = Z[[k]], id = id[[k]])
		
			## Update random effects conditional on fixed
			#do.inv.D <- try(solve(D[[k]] + diag(x=1e-5, nrow = num.ran[k])), silent=TRUE); 
			do.inv.D <- ginv(D[[k]]) # ## Need this when rows of D are zero
			
			if(pen.type[2] %in% c("adl","lasso")) 
				D2 <- diag(x=length(y)*lambda[2]/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num.ran[k]) 
			if(pen.type[2] == "scad") 
				D2 <- diag(x=length(y)*scad.deriv(sqrt(colSums(b[[k]]^2)),lambda=lambda[2],a=control$scad.a)/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num.ran[k]) 
			if(pen.type[2] == "mcp") 
				D2 <- diag(x=length(y)*mcp.deriv(sqrt(colSums(b[[k]]^2)),lambda=lambda[2],gamma=control$mcp.gamma)/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num.ran[k]) 
			if(!is.null(pen.weights$ran)) diag(D2) <- diag(D2)*c(pen.weights$ran[[k]])
			D2 <- do.inv.D + D2
			D2[which(D2 > 1e8)] <- 1e8
				
				
			if(family$family[1] == "gaussian") 
				new.bvec <- solve(t(bigZ)%*%bigZ + kronecker(diag(x=length.uniids[k]),phi*D2))%*%t(bigZ)%*%c(y-X%*%new.beta-all.Zb-offset) 
			if(family$family[1] == "LOGNO") 
				new.bvec <- solve(t(bigZ)%*%bigZ + kronecker(diag(x=length.uniids[k]),phi*D2))%*%t(bigZ)%*%c(log(y)-X%*%new.beta-all.Zb-offset) 
			if(!(family$family[1] %in% c("gaussian","LOGNO"))) {
				do.bvec <- try(solve(t(bigZ)%*%diag(x=new.weights)%*%bigZ + kronecker(diag(x=length.uniids[k]),D2))%*%t(bigZ)%*%diag(x=new.weights)%*%c(new.workresp-X%*%new.beta-all.Zb-offset), silent=TRUE) 

				if(!inherits(do.bvec,"try-error")) new.bvec <- do.bvec
				if(inherits(do.bvec,"try-error")) 
					new.bvec <- ginv(t(bigZ)%*%diag(x=new.weights)%*%bigZ + kronecker(diag(x=length.uniids[k]),D2))%*%t(bigZ)%*%diag(x=new.weights)%*%c(new.workresp-X%*%new.beta-all.Zb-offset) 
				}
			new.b[[k]] <- matrix(new.bvec, nrow = length.uniids[k], ncol = num.ran[k], byrow = TRUE)		
			sel.nonzero.b <- which(sqrt(colSums(new.b[[k]]^2)) > length.uniids[k]*conv.eps)
			if(length(sel.nonzero.b) > 0) new.b[[k]][,-sel.nonzero.b] <- 0
			if(length(sel.nonzero.b) == 0) new.b[[k]][,] <- 0
			
			
			## Update random effects covariance matrix 
			## Iterative estimator based on maximizing the Laplace approximation 
			## If multiple Z are applied, then the Laplace approximated estimator is done conditonally on all other random effects
			new.D[[k]] <- matrix(0,num.ran[k],num.ran[k])
			if(length(sel.nonzero.b) == 0) next;

			if(family$family[1] %in% c("gaussian","LOGNO")) 
				new.D[[k]][sel.nonzero.b,sel.nonzero.b] <- solve(t(as.matrix(Z[[k]][,sel.nonzero.b]))%*%as.matrix(Z[[k]][,sel.nonzero.b])/phi + do.inv.D[sel.nonzero.b,sel.nonzero.b]) + cov.wt(as.matrix(new.b[[k]][,sel.nonzero.b]), center = FALSE)$cov*(length.uniids[k]-1)/length.uniids[k]

			if(!(family$family[1] %in% c("gaussian","LOGNO"))) {
				for(i in 1:length.uniids[k]) { 
					sel.i <- which(id[[k]] == i)
 					new.D[[k]][sel.nonzero.b,sel.nonzero.b] <- new.D[[k]][sel.nonzero.b,sel.nonzero.b] + solve(t(as.matrix(Z[[k]][sel.i,sel.nonzero.b]))%*%diag(x=new.weights[sel.i])%*%as.matrix(Z[[k]][sel.i,sel.nonzero.b]) + do.inv.D[sel.nonzero.b,sel.nonzero.b]) + as.matrix(new.b[[k]][i,sel.nonzero.b])%*%t(as.matrix(new.b[[k]][i,sel.nonzero.b])) 
					}
				new.D[[k]] <- new.D[[k]]/length.uniids[k]
				}			

# 			if(ran.cov.est == "bb") { ## Based on (1/n_k) sum_{i=1}^{n_k} b_{ik}%*%t(b_{ik})
# 				new.D[[k]] <- cov.wt(new.b[[k]], center = FALSE, method = "ML")$cov
# 				}
# 				
# 			if(ran.cov.est == "wbb") { ## Based on (1/n_k) sum_{i=1}^{n_k} w_{ik} b_{ik}%*%t(b_{ik}), where w_{ik} is based on cluster size
# 				clus.size <- unlist(as.vector(table(id[[k]])))
# 				new.D[[k]] <- cov.wt(new.b[[k]], center = FALSE, method = "ML", wt = clus.size)$cov
# 				}
				
			if(length(sel.nonzero.b) == 0) { ranef.fullshrunk[k] <- 1 }
			rm(bigZ, new.bvec)
			}	

			
		## Solve scale parameters if required
		all.Zb <- numeric(length(y)) 
		for(k2 in 1:length(id)) { all.Zb <- all.Zb + as.numeric(make.bigZ(Z = Z[[k2]], id = id[[k2]])%*%c(t(new.b[[k2]]))) }
		if(family$family[1] == "gaussian") new.phi <- mean((y-X%*%new.beta-all.Zb-offset)^2)
		if(family$family[1] == "LOGNO") new.phi <- mean((log(y)-X%*%new.beta-all.Zb-offset)^2)

		if(family$family[1] == "negative.binomial") {
			new.phi <- try(1/theta.ml(y = y, mu = family$linkinv(X%*%new.beta+all.Zb+offset), limit = 100), silent = TRUE)
			if(inherits(new.phi,"try-error")) { new.phi <- phi }
			}
			
		if(family$family[1] == "Gamma") {
			prof.logl <- function(a, y, mu) {
				out <- -y*a/mu + (a-1)*log(y) - lgamma(a) - a*log(mu) + a*log(a)
				return(sum(out)) 
				}
			update.shape <- try(optimize(f = prof.logl, interval = c(1e-3,1e3), y = y, mu = family$linkinv(X%*%new.beta+all.Zb+offset), maximum = TRUE), silent = TRUE)
			if(inherits(update.shape,"try-error")) { new.shape <- shape }			
			if(!inherits(update.shape,"try-error")) { new.shape <- update.shape$max }
			}
			
		if(family$family[1] == "ZIP") {
			prof.logl <- function(a, y, mu) {
				out <- sum(dZIP(x=y, mu=mu, sigma=a, log=TRUE)) 
				return(out)
				}
			update.zeroprob <- try(suppressWarnings(optimize(f = prof.logl, interval = c(1e-4,0.9999), y = y, mu = family$mu.linkinv(X%*%new.beta+all.Zb+offset), maximum = TRUE)), silent = TRUE)
			if(inherits(update.zeroprob,"try-error")) { new.zeroprob <- zeroprob }
			if(!inherits(update.zeroprob,"try-error")) { new.zeroprob <- update.zeroprob$max }
			}
			

			
		diff <- sum((new.beta-beta)^2) 
		for(k in 1:length(id)) { diff <- diff + sum((new.D[[k]][upper.tri(new.D,diag=TRUE)]-D[[k]][upper.tri(D,diag=TRUE)])^2) }
		#for(k in 1:length(id)) { diff <- diff + sum((new.b[[k]]-b[[k]])^2) }
		if(control$trace) cat("Iteration:", counter,"   Error:", round(diff,4), "\n")
		if(control$trace) { 
			new.beta2 <- new.beta; names(new.beta2) <- colnames(X)
			print("Fixed effects:"); print(c(new.beta2)); #print(c(newbeta.l2))
			print("Covariance Matrices:"); print(new.D); 
			rm(new.beta2) 
			} 

			
		beta <- new.beta
		b <- new.b; 
		D <- new.D; 
		shape <- new.shape
		zeroprob <- new.zeroprob
		phi <- new.phi
		counter <- counter + 1;
		true.counter <- true.counter + 1; 
		if(counter < 5) true.counter <- 0

		
 		if(any(unlist(lapply(D, function(x) any(abs(x) > 30)) == TRUE)) || sum(abs(new.beta)) < 1e-3) {
 			if(control$trace) cat("Convergence issues. Restarting...\n")
 			restart.counter <- restart.counter + 1
			if(init.fit.run == 0) {
				get.init <- start.fixed(y = y, X = X, Z = Z, id = id, family = family, offset = offset, trial.size = trial.size)
				start <- build.start.fit.rpql(fit = get.init, id = id, num.ran = num.ran, cov.groups = cov.groups)
				init.fit.run <- 1
				}
 			
			new.beta <- beta <- start$fixef
			new.b <- b <- vector("list",length(id))
			for(k in 1:length(b)) { 
				new.b[[k]] <- b[[k]] <- matrix(rnorm(length.uniids[k]*num.ran[k], mean=0, sd = 0.2+0.2*(family$family!="poisson")),length.uniids[k],num.ran[k]) 
				}
			new.D <- D <- vector("list",length(id))
			for(k in 1:length(D)) { new.D[[k]] <- D[[k]] <- cov(new.b[[k]]) }

			new.phi <- phi <- new.shape <- shape <- new.zeroprob <- zeroprob <- 1
			if(family$family[1] == "gaussian") new.phi <- phi <- mean((y-X%*%beta-offset)^2) 
			if(family$family[1] == "LOGNO") new.phi <- phi <- mean((log(y)-X%*%beta-offset)^2) 
			if(family$family[1] == "negative.binomial") new.phi <- phi <- 0.01 
			if(family$family[1] == "Gamma") new.shape <- shape <- 1
			if(family$family[1] == "ZIP") new.zeroprob <- zeroprob <- 0.1
			
 			diff <- 1e4; cw.rpql.loglik <- 0; counter <- 0 
 			}		

		if(restart.counter > control$restarts) { 
			if(control$trace) stop("Sorry, but convergence could not be reached within the allocated number of control$restarts. Consider increasing lambda[2].\n") 
			}

		}
	## DONE!

		
	sel.nonzero.b <- sel.zero.b <- vector("list",length(id))
	for(k in 1:length(id)) { 
		sel.nonzero.b[[k]] <- which(sqrt(colSums(new.b[[k]]^2)) > length.uniids[k]*conv.eps) 
		if(length(sel.nonzero.b[[k]]) > 0) new.b[[k]][,-sel.nonzero.b[[k]]] <- 0
		}
	
	
	## The PQL likelihood
	new.etas <- X%*%new.beta + offset
	for(k2 in 1:length(id)) { new.etas <- new.etas + as.numeric(make.bigZ(Z = Z[[k2]], id = id[[k2]])%*%c(t(new.b[[k2]]))) }

	if(family$family[1] == "gaussian") { pql.loglik <- sum(dnorm(y, mean=new.etas, sd=sqrt(new.phi), log=TRUE)) } 
	if(family$family[1] == "Gamma") { pql.loglik <- sum(dgamma(y, shape=new.shape, scale=family$linkinv(new.etas)/new.shape, log=TRUE)) } 
	if(family$family[1] == "binomial") { pql.loglik <- sum(dbinom(y, size=trial.size, prob=family$linkinv(new.etas), log=TRUE)) } 
	if(family$family[1] == "poisson") { pql.loglik <- sum(dpois(y, lambda=family$linkinv(new.etas), log=TRUE)) } 
	if(family$family[1] == "negative.binomial") { pql.loglik <- sum(dnbinom(y, mu=family$linkinv(new.etas), size=1/new.phi, log=TRUE)) }
	if(family$family[1] == "LOGNO") { pql.loglik <- sum(dlnorm(y, meanlog=new.etas, sdlog=sqrt(new.phi), log=TRUE)) }
	if(family$family[1] == "ZIP") { pql.loglik <- sum(dZIP(y, mu=family$mu.linkinv(new.etas), sigma=new.zeroprob, log=TRUE)) }
	loglik1 <- pql.loglik
	
	for(k in 1:length(id)) { if(length(sel.nonzero.b[[k]]) > 0) 
		pql.loglik <- pql.loglik - 0.5*sum(rowSums(new.b[[k]][,sel.nonzero.b[[k]]]*(new.b[[k]][,sel.nonzero.b[[k]]]%*%solve(new.D[[k]][sel.nonzero.b[[k]],sel.nonzero.b[[k]]])))) 
		}
	names(sel.nonzero.b) <- names(id)
		
		
	## Information Criterion...argh!
	aic.v1 <- -2*loglik1 + 2*(sum(new.beta != 0) + sum(sapply(new.b, function(x) sum(x!=0))))
	bic.v2 <- -2*loglik1 + log(length(y))*sum(new.beta!=0) + sum(log(length.uniids)*sapply(new.D, function(x) sum(x[upper.tri(x,diag=T)]!=0))) ## Number of random effects to penalize is just number of non-zero elements in covariance matrix
	bic.v3 <- -2*loglik1 + log(length(y))*sum(new.beta!=0) + sum(log(length(y))*sapply(new.D, function(x) sum(x[upper.tri(x,diag=T)]!=0))) ## Number of random effects to penalize is just number of non-zero elements in covariance matrix
	hic.v1 <- -2*loglik1 + log(length(y))*sum(new.beta!=0) + 2*sum(sapply(new.b, function(x) sum(x!=0)))
	hic.v2 <- -2*loglik1 + log(length(y))*sum(new.beta!=0) + sum(sapply(new.b, function(x) sum(x!=0)))
	
	
	dof <- sum(new.beta != 0) + sum(sapply(new.b, function(x) sum(x!=0)))
	ics <- c(dof, aic.v1, bic.v2, bic.v3, hic.v1, hic.v2); 
	names(ics) = c(
		"# of estimated parameters",
		"AIC: 2*sum(beta!=0) + 2*sum(b!=0)",
		"BIC: log(nrow(X))*sum(beta!=0) + log(number of clusters)*sum(D!=0)",
		"BIC: log(nrow(X))*sum(beta!=0) + log(nrow(X))*sum(D!=0)",
		"Hybrid IC: log(nrow(X))*sum(beta!=0) + 2*sum(b!=0)",
		"Hybrid IC: log(nrow(X))*sum(beta!=0) + sum(b!=0)")

		
	## Garnishing
	out.list <- list(fixef = as.vector(new.beta), ranef = new.b, ran.cov = new.D, pql.logLik = pql.loglik, logLik1 = loglik1, phi = new.phi, shape = new.shape, zeroprob = new.zeroprob, family = family, n = length(y), trial.size = trial.size, id = id, lambda = lambda, pen.weights = pen.weights, pen.type = pen.type, ics = ics, nonzero.fixef = which(as.vector(new.beta)!=0), nonzero.ranef = sel.nonzero.b)

	if(save.data) { 
		out.list$y <- y; out.list$X <- X; outlist$Z <- Z; out.list$offset <- offset
		print("Please note the column names in Z have been overwritten for easier reference in the estimation algorithm. Apologies in advance") 
		}

		
	if(hybrid.est) {
		if(family$family[1] %in% c("ZIP","LOGNO")) { 
			print("Hybrid estimation not possible as lme4 does not do this. Sorry! Maybe try glmmADMB or MCMCglmm?")
			out.list$hybrid <- NULL 
			}
			
		else {	
			if(control$trace) 
				print("Performing hybrid estimation with lme4...might take a bit of time, and the labels for the parameters will be wrong. Apologies in advance!")

			make.dat <- data.frame(y, do.call(cbind,Z), X, do.call(cbind,id))
			restring1 <- NULL
			if(length(out.list$nonzero.fixef) > 0) restring1 <- paste(paste(colnames(X)[out.list$nonzero.fixef],collapse="+"),"-1")
			for(k in 1:length(id)) { 
				if(length(out.list$nonzero.ranef[[k]])) restring1 <- paste(restring1, "+ (",paste0(colnames(Z[[k]])[out.list$nonzero.ranef[[k]]],collapse="+"),"-1|",names(id)[k],")")
				}
			restring1 <- reformulate(restring1, response = "y")
	# 		print(restring1)
			if(family$family[1] == "gaussian") 
				final.fit <- suppressWarnings(lmer(restring1, REML = TRUE, offset = offset, data = make.dat))
			if(family$family[1] != "gaussian" & family$family[1] != "negative.binomial") 
				final.fit <- suppressWarnings(glmer(restring1, family = family, offset = offset, data = make.dat))
			if(family$family[1] == "negative.binomial") 
				final.fit <- suppressWarnings(glmer.nb(restring1, family = family, offset = offset, data = make.dat))

			out.list$hybrid <- final.fit ## Access covariance matrix as VarCorr(final.fit)[[1]][1:nrow(VarCorr(final.fit)[[1]]),1:nrow(VarCorr(final.fit)[[1]])]
			
			rm(make.dat)
			}
		}	
	

	names(out.list$fixef) <- colnames(X)
	names(out.list$ranef) <- names(out.list$ran.cov) <- names(id)
	for(k in 1:length(id)) { 
		rownames(out.list$ranef[[k]]) <- 1:length.uniids[k]; colnames(out.list$ranef[[k]]) <- colnames(Z[[k]])
		rownames(out.list$ran.cov[[k]]) <- colnames(out.list$ran.cov[[k]]) <- colnames(Z[[k]]) 
		}
		
	class(out.list) <- "rpql"
	out.list$call <- match.call()
	
	return(out.list)
	}

	
	
print.rpql <- function(x, ...) {
 	message("Call:")
 	print(x$call)
 	message()

 	cat("Family:", x$family$family[1], "\nPenalty type:", x$pen.type, "\nTuning parameters:", as.numeric(x$lambda), "\n") 
 	cat("Total number of observations:", x$n, "\n IDs (groups):", paste(names(x$id),sapply(x$id,function(a) length(unique(a))), sep = ": ", collapse = ";"), "\n") 
 	cat("Non-zero fixed effects:", x$nonzero.fixef, "\nNon-zero random effects by IDs (groups):", paste(names(x$id),x$nonzero.ranef, sep = ": ", collapse = ";"), "\n") 
 	#message("Information Criterion:")
 	#print(x$ics) 
 	}
	
	
print.summary.rpql <- function(x, ...) {
 	message("Call:")
 	print(x$call)
 	message()
 	
 	cat("Family:", x$family$family[1], "\nPenalty type:", x$pen.type, "\n Tuning parameters:", as.numeric(x$lambda), "\n") 
 	cat("Value of PQL at convergence:", x$logLik, "\n\n")
	#cat("Estimated model\n\t Non-zero fixed effects --", x$nonzero.fixef, "\n\t Non-zero random effects --", paste(names(x$id),x$nonzero.ranef, sep = ": ", collapse = "; "), "\n\n") 
  	message("Estimates of fixed effects:")
  	print(x$fixef); 
  	message()

  	message("Estimates of variances for random effects (diagonal elements of the covariance matrices):")
 	for(k in 1:length(x$ranef)) { message(names(x$id)[k], ":"); print(diag(x$ran.cov[[k]])); message() } 

 	if(x$family$family[1] %in% c("gaussian","lognormal")) message("Variance par ameter: ", x$phi, "\n")
 	if(x$family$family[1] %in% c("negative.binomial")) message("overdispersion parameter (V = mu + phi*mu^2): ", x$phi, "\n")
 	if((x$family$family[1] == "Gamma")) message("Shape parameter (V = mu^2/shape): ", x$shape, "\n")
 	if((x$family$family[1] == "ZIP")) message("Probability of structural zero: ", x$zeroprob, "\n") 
	}	
	
	
summary.rpql <- function(object, ...) {
 	gather.output <- list(call = object$call, fixef = round(object$fixef,3), ranef = lapply(object$ranef,round,3), ran.cov = lapply(object$ran.cov,round,3), logLik = round(object$logLik,3), family = object$family)
 	
 	if(object$family$family[1] %in% c("gaussian","lognormal","negative.binomial")) gather.output$phi <- round(object$phi,3)
 	if((object$family$family[1] == "Gamma")) gather.output$shape <- round(object$shape,3)
 	if((object$family$family[1] == "ZIP")) gather.output$zeroprob <- round(object$zeroprob,3)
 
 	gather.output$pen.type <- object$pen.type
 	gather.output$lambda <- object$lambda
 	gather.output$ics <- object$ics
	gather.output$id <- object$id 
 	gather.output$nonzero.fixef <- object$nonzero.fixef
 	gather.output$nonzero.ranef <- object$nonzero.ranef
 	gather.output$ics <- object$ics 
 	
 	class(gather.output) <- "summary.rpql"
 	gather.output 
 	}	
	
	
	