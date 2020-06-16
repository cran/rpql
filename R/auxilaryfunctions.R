#################
## Auxilary functions
##################
build.start.fit <- function(lme4.fit, id = NULL, gamma = 0, cov.groups = NULL) {
	if(!(class(lme4.fit)[1] %in% c("glmerMod","lmerMod"))) 
          stop("lme4.fit must be fit from the lme4 package. Thanks")

	out <- list(fixef = fixef(lme4.fit)+1e-4)
	get_ranefs <- as.list(ranef(lme4.fit))
	for(k in 1:length(get_ranefs)) 
          get_ranefs[[k]] <- as.matrix(get_ranefs[[k]])
	if(!is.null(id)) 
          names(get_ranefs) <- names(id)
	out$ranef <- get_ranefs
	out$cov.groups <- cov.groups
	
	get_varcors <- vector("list", length(get_ranefs))
     names(get_varcors) <- names(get_ranefs)
	for(k in 1:length(get_ranefs)) 
		get_varcors[[k]] <- VarCorr(lme4.fit)[[k]]+1e-4
	out$ran.cov <- get_varcors

	
	if(length(gamma) == 1) 
		gamma <- rep(gamma,2)
	#cat("Building adaptive lasso weights...\n")
	if(is.null(cov.groups)) 
          out$pen.weights <- list(fixed = 1/abs(out$fixef)^gamma[1])
	if(!is.null(cov.groups)) {
		make_adlweights <- (1/unlist(lapply(split(out$fixef,cov.groups), l2.norm)))[cov.groups] ## Split by cov.groups, calc L2 norm for each, convert to vector, then expand back...awesome!
		names(make_adlweights) <- names(out$fixef)
		out$pen.weights <- list(fixed = make_adlweights^gamma[1]) 
		}

	out$pen.weights$random <- vector("list", length(get_ranefs))
	names(out$pen.weights$random) <- names(get_ranefs)
	for(k in 1:length(get_ranefs)) 
          out$pen.weights$random[[k]] <- 1/diag(out$ran.cov[[k]]^gamma[2])
	
	return(out)
	}


## The marginal logL is calculated using Monte Carlo integration, exploiting the fact that conditional on all the random effects, all obs are independent. So generate a whole bunch of b's, then for b calculate the likelihood, then average over the b's
calc.marglogL <- function(new.data, fit, B = 1000) {
	if(!all(c("y","X","Z") %in% attributes(new.data)$names)) 
		stop("new.data must be a list containing y, X, Z")
	new.data$X <- as.matrix(new.data$X)
	n <- nrow(new.data$X)
	family <- fit$family

	if(any(sapply(new.data$Z,nrow) != n)) 
		stop("The number of rows in each element of new.data$Z should be equal to the number of rows in new.data$X. Thanks.")

	get_ncolZ <- sapply(new.data$Z,ncol)
	get_ncolranef_covmat <- sapply(fit$ran.cov,ncol)
	if(any(get_ncolZ != get_ncolranef_covmat)) 
		stop("The number of columns in each in the list new.data$Z should equal to dimension of the corresponding element in the list fit$ran.cov. Thanks.")
	rm(get_ncolZ, get_ncolranef_covmat)
	for(k in 1:length(fit$ran.cov)) 
          new.data$Z[[k]] <- as.matrix(new.data$Z[[k]])

	
	loglik <- numeric(B); 
	## Generate random effects
	get_nonzeroranef_covmat <- genb <- vector("list",length(fit$ran.cov)); 
	for(k in 1:length(fit$ran.cov)) { 
		get_nonzeroranef_covmat[[k]] <- which(diag(fit$ran.cov[[k]]) > 0)
		genb[[k]] <- matrix(0, B, ncol(new.data$Z[[k]]))
		genb[[k]][,get_nonzeroranef_covmat[[k]]] <- rmvnorm(B, mean = rep(0,length(get_nonzeroranef_covmat[[k]])), as.matrix(fit$ran.cov[[k]][get_nonzeroranef_covmat[[k]],get_nonzeroranef_covmat[[k]]]))
		}


	eta <- new.data$X %*% fit$fixef + fit$offset
	for(k0 in 1:B) {
		eta2 <- eta
		for(k in 1:length(fit$ran.cov)) 
               eta2 <- eta2 + new.data$Z[[k]] %*% genb[[k]][k0,]
		
		if(family$family[1] == "gaussian") 
               tmp_loglik <- dnorm(new.data$y, mean = eta2, sd=sqrt(fit$phi), log = TRUE)
		if(family$family[1] == "Gamma") 
               tmp_loglik <- (dgamma(new.data$y, shape = fit$shape, scale = family$linkinv(eta2)/fit$shape, log = TRUE))
		if(family$family[1] == "poisson") 
               tmp_loglik <- (dpois(new.data$y, lambda=family$linkinv(eta2), log = TRUE))
		if(family$family[1] == "binomial") 
               tmp_loglik <- (dbinom(new.data$y, size=fit$trial.size, prob=family$linkinv(eta2), log = TRUE))
		if(family$family[1] == "negative.binomial") 
               tmp_loglik <- (dnbinom(new.data$y, mu=family$linkinv(eta2), size=1/fit$phi, log = TRUE))
		if(family$family[1] == "LOGNO") 
               tmp_loglik <- (dlnorm(new.data$y, meanlog=eta, sdlog=sqrt(fit$phi), log=TRUE))
		if(family$family[1] == "ZIP") 
               tmp_loglik <- (dZIP(new.data$y, mu=family$mu.linkinv(eta2), sigma=fit$zeroprob, log=TRUE))

		tmp_loglik[!is.finite(tmp_loglik)] <- NA
		loglik[k0] <- exp(sum(tmp_loglik,na.rm=TRUE))
		} 
	
	return(mean(loglik[is.finite(loglik)]))
	}


	
## Dataset generation for GLMM
## Intercept must be manually included in X and Z if desired
#id = list(cluster = rep(1:n,each=m), cluster2 = id2); beta = truebeta; D = list(cluster = true.D, cluster2 = true.D2); trial.size = 1; family = "binomial"; phi = NULL; upper.count = Inf
gendat.glmm <- function(id, X, beta, Z, D, trial.size = 1, family = gaussian(), phi = NULL, shape = NULL, zeroprob = NULL, upper.count = Inf) {
	X <- as.matrix(X)
	n <- nrow(X)
	if(ncol(X) != length(beta)) 
		stop("The number of columns in X should be equal to the length of beta. Thanks")
	
	if(sum(X[,1]) != nrow(X)) 
		warning("Has an intercept column being included in X?")
	if(!is.list(id) || !is.list(D) || !is.list(Z)) 
		stop("Please supply id, Z, D, as lists. The elements in the three lists corresponding to a vector of IDs, and the model matrix and the random effects covariance matrix for those IDs, respectively. Please see the help file for further details. ")
	if(length(unique(length(id),length(D),length(Z))) != 1) 
		stop("The number of elements in the lists id, D, and Z should be the same. Thanks.")

	get_dimD <- sapply(D, ncol)
	get_ncolZ <- sapply(Z, ncol)
	if(any(get_dimD!=get_ncolZ)) 
		stop("The dimension of each element in the list D should be equal to the number of columns in the corresponding element of list Z. Thanks.")

		
	get_lengths_id <- sapply(id,length)
	if(any(get_lengths_id!=nrow(X))) 
		stop("The length of each element in the list id should be equal to the number of rows in X and Z. Thanks.")
	rm(get_lengths_id, get_dimD, get_ncolZ)
 
	
	if(!(family$family[1] %in% c("gaussian","poisson","binomial","negative.binomial","Gamma","LOGNO","ZIP"))) 
		stop("Current version does not permit specified family. Sorry!")
	if(family$family[1] == "Gamma" & is.null(shape)) 
		stop("Please supply shape parameter for the Gamma family (Variance = mu^2/shape).") 
	if(family$family[1] == "negative.binomial" & is.null(phi)) 
		stop("Please supply the overdispersion parameter phi for the negative binomial family; note the variance is parameterized as V = mu + phi*mu^2.") 
	if(family$family[1] == "gaussian" & is.null(phi)) 
		stop("Please supply the variance parameter phi for the Gaussian family; note the variance is parameterized as V = phi.") 
	if(family$family[1] == "LOGNO" & is.null(phi)) 
		stop("Please supply the variance parameter phi for the lognormal family; note the variance is parameterized as V = phi on the log scale).") 
	if(family$family[1] == "ZIP" & is.null(zeroprob)) 
		stop("Please supply zeroprob, the probability for obtaining a structural zero, for the ZIP family.") 
	
	if(is.null(colnames(X))) 
          colnames(X) <- paste("x",1:ncol(X),sep="")
	for(k in 1:length(Z)) 
          colnames(Z[[k]]) <- paste("z",k,1:ncol(Z[[k]]),sep="")

	get.n.id <- numeric(length(id))
	for(k in 1:length(id)) 
          get.n.id[k] <- length(unique(id[[k]]))
	
	for(k in 1:length(D)) { 
          D[[k]] <- as.matrix(D[[k]])
          Z[[k]] <- as.matrix(Z[[k]]) 
          }
	for(k in 1:length(id)) 
         { 
          id[[k]] <- as.integer(id[[k]]) 
          }

	
	## Generate random effects
	get_nonzeroD <- trueb <- vector("list",length(D)); 
	for(k in 1:length(trueb)) { 
		get_nonzeroD[[k]] <- which(diag(D[[k]]) > 0)
		trueb[[k]] <- matrix(0, nrow = get.n.id[k], ncol = ncol(Z[[k]]))
		trueb[[k]][,get_nonzeroD[[k]]] <- rmvnorm(get.n.id[k], rep(0,length(get_nonzeroD[[k]])), as.matrix(D[[k]][get_nonzeroD[[k]],get_nonzeroD[[k]]]))
		}

		
	## Generate response	
	sim_y <- numeric(n)
	eta <- X %*% beta 
	for(k in 1:length(trueb)) 
          eta <- eta + rowSums(Z[[k]]*trueb[[k]][id[[k]],])

	if(family$family[1] == "gaussian") 
          sim_y <- rnorm(n, mean = eta, sd=sqrt(phi))
	if(family$family[1] == "LOGNO") 
          sim_y <- rlnorm(n, meanlog=eta, sdlog=sqrt(phi))
	if(family$family[1] == "Gamma") 
          sim_y <- rgamma(n, shape=shape, scale=family$linkinv(eta)/shape)
	if(family$family[1] == "binomial") 
          sim_y <- rbinom(n, size=trial.size, prob=family$linkinv(eta))
	for(i in 1:n) {
		if(family$family[1] %in% c("poisson","negative.binomial","ZIP")) { 
			if(family$family[1] == "ZIP") 
                    testsim_y <- rZIP(1, mu=family$mu.linkinv(eta[i]), sigma=zeroprob)
			if(family$family[1] == "poisson") 
                    testsim_y <- rpois(1, lambda=family$linkinv(eta[i]))
			if(family$family[1] == "negative.binomial") 
                    testsim_y <- rnbinom(1, mu=family$linkinv(eta[i]), size=1/phi)
			try_counter <- 0
			
			while(testsim_y > upper.count & try_counter < 5000) {
				if(family$family[1] == "ZIP") 
                         testsim_y <- rZIP(1, mu=family$mu.linkinv(eta[i]), sigma=zeroprob)
				if(family$family[1] == "poisson") 
                         testsim_y <- rpois(1, lambda=family$linkinv(eta[i]))
				if(family$family[1] == "negative.binomial") 
                         testsim_y <- rnbinom(1, mu=family$linkinv(eta[i]), size=1/phi)
				try_counter <- try_counter + 1
				}
			sim_y[i] <- testsim_y
			}		
		}

		
	out <- list(y = sim_y, id = id, X = X, Z = Z, beta = beta, b = trueb, D = D, phi = phi, shape = shape, zeroprob = zeroprob, trial.size = trial.size, family = family, nonzero.beta = which(beta!=0), nonzero.b = get_nonzeroD)
	names(out$b) <- names(out$nonzero.b) <- names(D)
	
	return(out)
	}

	
	
l2.norm <- function(x) 
     sqrt(sum(x^2)) 


lseq <- function (from, to, length, decreasing = FALSE) {
	stopifnot(from > 0)
	out <- 10^(seq(log10(from), log10(to), length.out = length))
	out <- out[order(out, decreasing = decreasing)]; 
	return(out)
	}

	
## First deriative of MC+ penalty
mcp.deriv <- function(x, lambda, gamma = 2) {
	x2 <- abs(x)
	help1 <- sapply(x2, function(theta) { lambda*(1 - theta/(gamma*(lambda+1e-8))) } )
	out <- help1*(x2 < gamma*lambda) + 0*(x2 >= gamma*lambda)
	return(out)
	}
	

## Create a negative binomial family
nb2 <- function() {
	link <- "log"
	linkfun <- function(mu) log(mu)
	linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
	mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
	variance <- function(mu,phi) mu+phi*mu^2
	valideta <- function(eta) TRUE
	validmu <- function(mu) all(mu > 0)
  
	structure(list(family = "negative.binomial", link = "log", linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance, valideta = valideta, validmu = validmu, name = link), class = "family")
	}


	    
## First deriative of SCAD penalty
scad.deriv <- function(x, lambda, a = 3.7) {
	x2 <- abs(x)
	help1 <- sapply(x2, function(theta) { max(a*lambda-theta,0)/((a-1)*(lambda+1e-8)) } )
	out <- lambda*((x2 <= lambda) + help1*(x2 > lambda))
	return(out)
	}


	
	
