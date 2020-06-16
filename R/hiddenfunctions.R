## Used only within rpql function
build.start.fit.rpql <- function(fit, id, num.ran, gamma = 0, cov.groups = NULL) {
	length_uniids <- numeric(length(id))
	for(k in 1:length(id)) { 
		length_uniids[k] <- length(unique(id[[k]])) 
		}
# 	print(class(fit))
		
	if(class(fit)[1] %in% c("glmerMod","lmerMod")) 
		get_fit <- build.start.fit(lme4.fit = fit, gamma = gamma, cov.groups = cov.groups)
	if(class(fit)[1] %in% c("glm","lm")) {
		get_fit <- list(fixef = fit$coefficients, ranef = vector("list",length(id)), ran.cov = vector("list",length(id)))
		for(k in 1:length(get_fit$ranef)) { 
			get_fit$ranef[[k]] <- matrix(0.01,length_uniids[k],num.ran[k]) 
			get_fit$ran.cov[[k]] <- diag(x = 1, nrow = num.ran[k]) 
			}
		names(get_fit$ranef) <- names(get_fit$ran.cov) <- names(id)	
		}
		
		
	return(get_fit)
	}

	
fillin.control <- function(control) {
	if(!("tol" %in% names(control))) 
          control$tol <- 1e-4
	if(!("maxit" %in% names(control))) 
          control$maxit <- 100
	if(!("trace" %in% names(control))) 
          control$trace <- FALSE
	if(!("restarts" %in% names(control))) 
          control$restarts <- 5
	if(!("scad.a" %in% names(control))) 
          control$scad.a <- 3.7
	if(!("mcp.gamma" %in% names(control))) 
          control$mcp.gamma <- 2
	if(!("seed" %in% names(control))) 
          control$seed <- NULL 
	if(!("lasso.lambda.scale" %in% names(control))) 
          control$lasso.lambda.scale <- TRUE

     return(control)
     }
	
	

## Generate starting values for the fixed effects by fitting a lme4 object with a small number of iterations
## For binomial and gaussian data, starting using a straight glm and setting ranef = list(cluster = matrix(0.01,n,9), ran.cov = list(cluster = diag(x = 0.5, nrow = 9))) works
start.fixed <- function(y, X, Z, id, family = gaussian(), offset = NULL, trial.size = 1) {
	make_dat <- data.frame(y, do.call(cbind,Z), X, do.call(cbind,id))
	restring1 <- paste(paste(colnames(X),collapse="+"),"-1")
	for(k in 1:length(id)) { 
		restring1 <- paste(restring1, "+ (",paste0(colnames(Z[[k]]),collapse="+"),"-1|",names(id)[k],")")
		}
	
	if(family$family[1] == "gaussian") {
		restring1 <- paste(paste(colnames(X),collapse="+"),"-1")
		restring2 <- reformulate(restring1, response = "y")		
		final_fit <- suppressWarnings(lm(restring2, offset = offset, data = make_dat))
		}
		
	if(family$family[1] == "binomial") {
		restring1 <- paste(paste(colnames(X),collapse="+"),"-1")
		restring2 <- as.formula(paste("cbind(y,trial.size-y) ~", restring1))
 		final_fit <- suppressWarnings(glm(restring2, family = family, offset = offset, data = make_dat))
# 		final_fit <- suppressWarnings(glmer(restring2, family = family, offset = offset, data = make_dat, control = glmerControl(optCtrl=list(maxfun=10))))
		}
	
	if(family$family[1] %in% c("ZIP","poisson")) {
		restring2 <- reformulate(restring1, response = "y")		
		final_fit <- suppressWarnings(glmer(restring2, family = poisson(), offset = offset, data = make_dat, control = glmerControl(optCtrl=list(maxfun=10))))
		}
		
	if(family$family[1] %in% c("Gamma")) {
		restring2 <- reformulate(restring1, response = "y")		
		final_fit <- suppressWarnings(glmer(restring2, family = Gamma(), offset = offset, data = make_dat, control = glmerControl(optCtrl=list(maxfun=10))))
		}

	if(family$family[1] %in% c("LOGNO")) {
		restring2 <- as.formula(paste("log(y) ~", restring1))
		final_fit <- suppressWarnings(lmer(restring2, REML = FALSE, offset = offset, data = make_dat, control = lmerControl(optCtrl=list(maxfun=10))))
		}
	
	if(family$family[1] %in% c("negative.binomial")) {
		restring2 <- reformulate(restring1, response = "y")		
		final_fit <- suppressWarnings(glmer.nb(restring2, offset = offset, data = make_dat, control = glmerControl(optCtrl=list(maxfun=10))))
		}

	return(final_fit)
	}



# ## Evaluate values of penalty functions
# pen.val <- function(beta,lambda,pen.type) {
# 	if(pen.type == "lasso") out <- sum(lambda*abs(beta))
# 	if(pen.type == "scad") { 
# 		a <- 3.7
# 		out <- sum(lambda*abs(beta)*(abs(beta)<lambda) - (beta^2 - 2*3.7*lambda*abs(beta) + lambda^2)/(2*(a-1))*(lambda < abs(beta) & abs(beta) < a*lambda) + lambda^2*(a+1)/2*(abs(beta) > lambda)) }
# 
# 	return(out)
# 	}
	
