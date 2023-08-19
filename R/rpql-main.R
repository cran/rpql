############
## Regularized PQL (with a penalty on fixed effects and a group penalty on the random effects
############

## Version 0.8.1
## - Fixed a minor bug on lasso.lambda.scale

## TODO: 1) consider using threshold operators, including LLA and threshold; 
## 2) Very slow for Gaussian responses!
## 3) Maybe scrap using C++ functions due to instability
##############
#  library(mvtnorm)
#  library(Matrix)
#  library(MASS)
#  library(lme4)
#  library(gamlss.dist)
#  library(Rcpp)
#  library(inline)
# 
# source("rpql08/R/hiddenfunctions.R")
# source("rpql08/R/auxilaryfunctions.R")
# sourceCpp(file = "rpql08/src/RcppExports.cpp")
# sourceCpp(file = "rpql08/src/newcoeffn.cpp")


## pen.weights is a list with up to two elements for adl; first contains weights for fixed, and second contains a list of weights for random. 
# y = simy$y; X = simy$X; Z = simy$Z; id = simy$id; family = binomial(); lambda = 0.1; pen.type = "mcp"; hybrid.est = FALSE; offset = NULL; trial.size = 1; pen.weights = NULL; cov.groups = NULL; trace = TRUE; control$restarts <- 10; control$seed = NULL; intercept = TRUE; save.data = FALSE; tol = 1e-4; control$maxit = 100; start = NULL

rpql <- function(y, ...) 
     UseMethod("rpql")


rpqlseq <- function(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, pen.type = "lasso", start = NULL, cov.groups = NULL, 
     pen.weights = NULL, offset = NULL, intercept = TRUE, save.data = FALSE, 
     control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, lasso.lambda.scale = TRUE, seed = NULL), ...) {
	
	lambda <- as.matrix(lambda)
	if(ncol(lambda) == 1) {
          lambda[,1] <- sort(lambda[,1], decreasing = FALSE)
          }
          
	
	control <- fillin.control(control) 
	if(!control$lasso.lambda.scale)
          warnings("Scaling factor for the second tuning parameter has been turned off for lasso and adaptive lassp penalties. 
               This may cause issues as the group coefficient penalty on the random effects is on a different scale to the single 
               coefficient penalty on the fixed effects.")
	
	## Do a starting fit for formatting reasons
	start_fit <- rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[1,], pen.type = "lasso", start = start) 
	collect_ics <- matrix(0, nrow = nrow(lambda), ncol = length(start_fit$ics)+1) ## First column in start_fit$ics is DoF
	colnames(collect_ics) <- c("PQL Likelihood", names(start_fit$ics)) 

	best_fits <- vector("list", length(start_fit$ics)-1); 
	names(best_fits) <- names(start_fit$ics)[-1]
	rm(start_fit)
			
	for(l1 in 1:nrow(lambda)) {
		message("Onto ", l1)
		if(l1 == 1) 
               prev_fit <- start
               
          new_fit <- rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[l1,], pen.type = pen.type, pen.weights = pen.weights, start = prev_fit, 
               cov.groups = cov.groups, hybrid.est = FALSE, offset = offset, intercept = intercept, save.data = save.data, control = control)
		
		if(l1 == 1){
            for(k in 1:length(best_fits)) 
                best_fits[[k]] <- new_fit
            }

		if(l1 > 1){
               for(l2 in 1:length(best_fits)) {
                    if(new_fit$ics[l2+1] < best_fits[[l2]]$ics[l2+1]) 
                         best_fits[[l2]] <- new_fit
                    }
               }
		prev_fit <- new_fit
		collect_ics[l1,] <- c(new_fit$logLik1, new_fit$ics)
		}
	
	out <- list(best.fits = best_fits, collect.ics = collect_ics, lambda = lambda, control = control)	
	}



#y = simy$y
#X = x
#Z = simy$Z
#id = simy$id
#family = binomial()
#lambda = 0
#pen.type = "lasso"
#trial.size = 10
#trace = TRUE
#start = gensatfit
#hybrid.est = FALSE; offset = NULL; 
#pen.weights = NULL; 
#cov.groups = NULL; 
#control = list(tol = 1e-4, maxit = 100, trace = TRUE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, seed = NULL); 
#intercept = TRUE; 
#save.data = FALSE; 

rpql.default <- function(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, pen.type = "lasso", start = NULL, cov.groups = NULL, 
     pen.weights = NULL, hybrid.est = FALSE, offset = NULL, intercept = TRUE, save.data = FALSE, 
     control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, lasso.lambda.scale = TRUE, seed = NULL), ...) {

     conv.eps <- 1e-3
     control <- fillin.control(control) 
	
     y <- as.vector(y)
     X <- as.matrix(X)
     num_fixed <- ncol(X)
     if(is.null(colnames(X))) 
          colnames(X) <- paste0("x",1:ncol(X)) 
     for(k in 1:length(Z)) 
          colnames(Z[[k]]) <- paste0("z",k,1:ncol(Z[[k]]))
     if(det(crossprod(X)) == 0) 
          warning("Is X rank-deficient? Please double check!", immediate. = TRUE) 
     if(!is.list(Z) || !is.list(id)) 
          stop("Please supply id and Z as lists. The element in the two lists correspond to a vector of IDs and the random effects model matrix for those IDs respectively.")
     if(sum(duplicated(names(id))) > 0 || sum(duplicated(names(Z))) > 0) 
          stop("Please ensure the names of the elements in the lists id (and Z) are unique.")
     if(length(Z) != length(id)) 
          stop("The number of elements in Z and id should be the same. Thanks.")
# 	if(any(apply(X,2,sd) != 1) || any(apply(X,2,sd) != 1))
# 		message("You may want to consider standardizing your X and Z matrix prior to penalization.")
	

          
     if(is.null(offset)) 
          offset <- rep(0,length(y))
     if(length(offset) != length(y)) 
          stop("Offset should have the same length as y. Thanks.")

     num_ran <- get_lengths_id <- get_nrows_Z <- length_uniids <- numeric(length(Z)); 
     for(k in 1:length(Z)) { 
          Z[[k]] <- as.matrix(Z[[k]])
          num_ran[k] <- ncol(Z[[k]])
          get_nrows_Z[k] <- nrow(Z[[k]]) 
          }
     for(k in 1:length(id)) { 
          id[[k]] <- as.integer(id[[k]])
          get_lengths_id[k] <- length(id[[k]])
          length_uniids[k] <- length(unique(id[[k]])) 
          }
     if((get_nrows_Z != get_lengths_id) || (get_nrows_Z != nrow(X)) || (nrow(X) != get_lengths_id)) 
          stop("The number of rows in X, the length of each element in list id, and the number of rows in each element in list Z should all be the same. Thanks.")
     rm(get_lengths_id, get_nrows_Z)

	
     if(!is.null(cov.groups)) {
          actual_cov_groups <- cov.groups 
          if(control$trace) 
               message("Group penalization performed as cov.groups has been supplied.")
          if((ncol(X) != length(cov.groups))) 
               stop("If supplied, the length of cov.groups should be equal to the number of columns in X. Thanks") 
          }
     if(is.null(cov.groups)) 
          actual_cov_groups <- 1:ncol(X)
               
               
     if(length(pen.type) == 1) 
          pen.type <- rep(pen.type,2)
     if(sum(pen.type %in% c("lasso","scad","adl","mcp")) != 2) 
          stop("Current version does not permit specified penalty. Sorry!")		
     if(pen.type[1] == "adl" & !is.vector(pen.weights$fixed)) 
          stop("Please supply weights for the fixed effects in adaptive lasso. Specifically, pen.weights$fixed should be a vector.")
     if(pen.type[2] == "adl" & !is.list(pen.weights$ran)) 
          stop("Please supply weights for the random effects in adaptive lasso. Specifically, pen.weights$ran should be a list with the same number of elements as lists Z and id.")
     if(!is.null(pen.weights$fixed) & length(pen.weights$fixed) != ncol(X)) 
          stop("Weights provided for fixed effects must have the same length as the number of columns in X. Thanks.")
     if(!is.null(pen.weights$ran)) 
          { 
          get_lengths_weights <- sapply(pen.weights$ran, length)
          if(!all(get_lengths_weights == num_ran)) 
               stop("The length of each element in pen.weights$ran should equal the number of columns in the corresponding element of list Z. Thanks.") 
          rm(get_lengths_weights) 
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
          
          
     if(length(lambda) == 1) {
          lambda <- rep(lambda, 2)
          if(!control$lasso.lambda.scale)
               warnings("Scaling factor for the second tuning parameter has been turned off for lasso and adaptive lassp penalties. 
                    This may cause issues as the group coefficient penalty on the random effects is on a different scale to the single 
                    coefficient penalty on the fixed effects.")
          }
     lambda <- as.vector(unlist(lambda))
     if(!(length(lambda) %in% c(1,2))) 
          stop("lambda should be a vector of length 1 or 2.")
     runif(1)
     old <- .Random.seed
     on.exit( { .Random.seed <<- old } )
     if(!is.null(control$seed)) 
          set.seed(control$seed)
               
     init_fit_run <- 0
     if(is.null(start)) {
          get_init <- start.fixed(y = y, X = X, Z = Z, id = id, family = family, offset = offset, trial.size = trial.size)
          start <- build.start.fit.rpql(fit = get_init, id = id, num.ran = num_ran, cov.groups = cov.groups)
          init_fit_run <- 1
          }

               
     new_beta <- beta <- start$fixef
     new_b <- b <- start$ranef
     if(!is.list(new_b)) 
          stop("start$ranef should be a list of random effects. Thanks.")
     for(k in 1:length(new_b)) 
          new_b[[k]] <- b[[k]] <- as.matrix(b[[k]])

     if(is.null(start$ran.cov)) { 
          new_ranef_covmat <- ranef_covmat <- vector("list",length(id))
          for(k in 1:length(ranef_covmat)) 
               new_ranef_covmat[[k]] <- ranef_covmat[[k]] <- cov(new_b[[k]]) + diag(x=rep(1e-4,num_ran[k])) 
          }
     if(!is.null(start$ran.cov)) 
          new_ranef_covmat <- ranef_covmat <- start$ran.cov
          

     ## Construcs big block diagonal matrix based on an a vector of id and a Z matrix. Basically converts Z to a block Z.
     make.bigZ <- function(Z, id) {
          matlist <- split(as.data.frame(Z), id)
          matlist <- lapply(matlist, function(x) as.matrix(x))
          bigZ <- matrix(0, nrow = nrow(Z), ncol = ncol(Z)*length(matlist))
          for(k in 1:length(matlist)) { 
               bigZ[which(id == k), (k*ncol(Z)-ncol(Z)+1):(k*ncol(Z))] <- matlist[[k]] 
               }
          bigZ <- Matrix(bigZ, sparse = TRUE)
          return(bigZ)
          }

     if(family$family[1] %in% c("LOGNO","gaussian")) 
          XtWX <- crossprod(X) ## Define these outside as they do not change with iteration
     all_bigZ <- vector("list", length(id))
     for(k in 1:length(id)) 
          all_bigZ[[k]] <- make.bigZ(Z = Z[[k]], id = id[[k]]) ## Define these outside as they do not change with iteration, but save them as sparse matrices

          
     new_phi <- phi <- new_shape <- shape <- new_zeroprob <- zeroprob <- 1
     all_Zb <- numeric(length(y)) 
     for(k2 in 1:length(id)) 
          all_Zb <- all_Zb + as.numeric(all_bigZ[[k2]] %*% c(t(b[[k2]])))
     if(family$family[1] == "gaussian") 
          new_phi <- phi <- mean((y - X %*% beta - all_Zb - offset)^2) 
     if(family$family[1] == "LOGNO") 
          new_phi <- phi <- mean((log(y) - X %*% beta - all_Zb - offset)^2) 
     if(family$family[1] == "negative.binomial") 
          new_phi <- phi <- 0.01 
     if(family$family[1] == "Gamma") 
          new_shape <- shape <- 1
     if(family$family[1] == "ZIP") 
          new_zeroprob <- zeroprob <- 0.1
     rm(all_Zb)

     if(family$family[1] %in% c("LOGNO","gaussian")) 
          XtWX <- crossprod(X) ## Define these outside as they do not change with iteration
     all_bigZ <- vector("list", length(id))
     for(k in 1:length(id)) 
          all_bigZ[[k]] <- make.bigZ(Z = Z[[k]], id = id[[k]]) ## Define these outside as they do not change with iteration, but save them as sparse matrices

          
     diff <- 1e4
     cw_rpqllogLik <- 0
     counter <- true.counter <- restart_counter <- 0
     ranef_fullshrunk <- rep(0,length(Z)) ## Indicator if random effects has been shrunk to zero. If 1, then ignore updating it!
    while(diff > control$tol & true.counter < control$maxit & restart_counter <= control$restarts) {
          all_Zb <- numeric(length(y)) 
          for(k2 in 1:length(id)) 
               all_Zb <- all_Zb + as.numeric(all_bigZ[[k2]] %*% c(t(b[[k2]])))

            
          if(!(family$family[1] %in% c("gaussian","LOGNO"))) {
               new_etas <- X %*% beta + all_Zb + offset
               new_etas[new_etas > 30] <- 30
               new_etas[new_etas < -30] <- -30

               if(family$family[1] == "negative.binomial") 
                    new_weights <- (family$mu.eta(new_etas))^2/(family$variance(family$linkinv(new_etas),phi=phi)) ## 1/(variance*g'(mu)^2)
               if(family$family[1] == "Gamma")
                    new_weights <- (family$mu.eta(new_etas))^2/(family$variance(family$linkinv(new_etas))/shape) 
               if(family$family[1] == "ZIP") {
                    new_weights <- (poisson()$mu.eta(new_etas))^2/(poisson()$variance(poisson()$linkinv(new_etas))) 
                    weightmulti <- (1-zeroprob)*dpois(y,lambda=family$mu.linkinv(new_etas))/dZIP(x=y,mu=family$mu.linkinv(new_etas),sigma=zeroprob)
                    weightmulti[!is.finite(weightmulti)] <- 0
                    new_weights <- new_weights*weightmulti
                    rm(weightmulti)
                    }
               if(!(family$family[1] %in% c("ZIP","negative.binomial","Gamma"))) 
                    new_weights <- trial.size*(family$mu.eta(new_etas))^2/(family$variance(family$linkinv(new_etas))) ## 1/(variance*g'(mu)^2)
               if(family$family[1] != "ZIP") 
                    new_workresp <- new_etas + (y/trial.size-family$linkinv(new_etas))/family$mu.eta(new_etas) ## z = eta + (y-mu)*g'(mu)
               if(family$family[1] == "ZIP") 
                    new_workresp <- new_etas + (y-family$mu.linkinv(new_etas))/poisson()$mu.eta(new_etas) 

               new_weights <- as.vector(new_weights)	
               new_weights[new_weights > 1e4] <- 1e4
               #XtWX <- crossprod(X*sqrt(new_weights)) ## REMOVE SINCE USING C++
               #XtWZ <- crossprod(X, (new_workresp-all_Zb-offset)*new_weights) ## REMOVE SINCE USING C++
               }

                         
          ## Update fixed effects conditional on random
          absbetas <- sapply(split(as.vector(beta),actual_cov_groups), l2.norm)[actual_cov_groups] + 1e-6
          penmat <- matrix(0,length(absbetas),length(absbetas))		
          if(pen.type[1] %in% c("adl","lasso")) 
               diag(penmat) <- length(y)*lambda[1]/as.vector(absbetas)
          if(pen.type[1] == "scad") 
               diag(penmat) <- length(y)*scad.deriv(absbetas, lambda=lambda[1], a = control$scad.a)/as.vector(absbetas)
          if(pen.type[1] == "mcp") 
               diag(penmat) <- length(y)*mcp.deriv(absbetas, lambda=lambda[1], gamma = control$mcp.gamma)/as.vector(absbetas)
          if(!is.null(pen.weights$fixed)) 
               diag(penmat) <- diag(penmat)*c(pen.weights$fixed)
          if(intercept) 
               diag(penmat)[1] <- 0 
          penmat[which(penmat > 1e8)] <- 1e8
        
          if(family$family[1] == "gaussian") 
               new_beta <- chol2inv(chol(XtWX + phi*penmat)) %*% crossprod(X, y-all_Zb-offset)
          if(family$family[1] == "LOGNO") 
                new_beta <- chol2inv(chol(XtWX + phi*penmat)) %*% crossprod(X, log(y)-all_Zb-offset)
          if(!(family$family[1] %in% c("gaussian","LOGNO"))) {
               new_beta <- newcoeffn(X*sqrt(new_weights), penmat, (new_workresp-all_Zb-offset)*sqrt(new_weights))			
            #new_beta <- try(chol2inv(chol(XtWX + penmat)) %*% XtWZ, silent=TRUE) ## REMOVE SINCE USING C++
#  			if(inherits(new_beta,"try-error")) 
#                      new_beta <- ginv(XtWX + penmat) %*% XtWZ
               }
          if(lambda[1] > 0) {
               if(intercept) 
                    new_betaint <- new_beta[1]
                    newbeta.l2 <- sapply(split(as.vector(new_beta),actual_cov_groups), l2.norm)
                    get.zeros <- which(newbeta.l2 < conv.eps)
                    new_beta[which(actual_cov_groups %in% get.zeros)] <- 0
               if(intercept) 
                    new_beta[1] <- new_betaint ## Do not allow intercept to be shrunk to zero
               rm(get.zeros, newbeta.l2)
               }
                    

          ## Update random effects conditional on fixed
          for(k in 1:length(id)) {
               if(ranef_fullshrunk[k] == 1) 
                    next;
                
               all_Zb <- numeric(length(y)); 
               if(length(Z) > 1) { 
                    for(k2 in (1:length(id))[-k]) 
                         all_Zb <- all_Zb + as.numeric(all_bigZ[[k2]] %*% c(t(b[[k2]]))) 
                    } ## Conditioning on all linear predictors for other ids, if is more than one id
               bigZ <- as.matrix(all_bigZ[[k]])
        
               do.inv.ranef_covmat <- try(solve(ranef_covmat[[k]] + diag(x=1e-8, nrow = num_ran[k])), silent=TRUE); 
               #do.inv.ranef_covmat <- ginv(ranef_covmat[[k]]) # ## Need this when rows of ranef_covmat are zero
                
               if(pen.type[2] %in% c("adl","lasso")) {
                    if(!control$lasso.lambda.scale)
                         penmat2 <- diag(x=length(y)*lambda[2]/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
                    if(control$lasso.lambda.scale)
                         penmat2 <- diag(x=length(y)*(lambda[2]/sqrt(length_uniids[k]))/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
                    }
               if(pen.type[2] == "scad") 
                    penmat2 <- diag(x=length(y)*scad.deriv(sqrt(colSums(b[[k]]^2)),lambda=lambda[2],a=control$scad.a)/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
               if(pen.type[2] == "mcp") 
                    penmat2 <- diag(x=length(y)*mcp.deriv(sqrt(colSums(b[[k]]^2)),lambda=lambda[2],gamma=control$mcp.gamma)/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
               if(!is.null(pen.weights$ran)) 
                    diag(penmat2) <- diag(penmat2)*c(pen.weights$ran[[k]])
               penmat2 <- do.inv.ranef_covmat + penmat2
               penmat2[which(penmat2 > 1e8)] <- 1e8
                        
                        
               if(family$family[1] == "gaussian") {
                    #new_bvec <- chol2inv(chol(crossprod(bigZ) + kronecker(diag(x=length_uniids[k]),phi*penmat2))) %*% crossprod(bigZ, c(y - X%*%new_beta - all_Zb - offset)) 
                    new_bvec <- newcoeffn(bigZ, kronecker(diag(x=length_uniids[k]),phi*penmat2), c(y - X %*% new_beta - all_Zb - offset))
                    }
               if(family$family[1] == "LOGNO") {
                    #new_bvec <- chol2inv(chol(crossprod(bigZ) + kronecker(diag(x=length_uniids[k]),phi*penmat2))) %*% crossprod(bigZ, c(log(y) - X%*%new_beta -all_Zb-offset))
                    new_bvec <- newcoeffn(bigZ, kronecker(diag(x=length_uniids[k]),phi*penmat2), c(log(y) - X %*% new_beta - all_Zb - offset))
                    }
               if(!(family$family[1] %in% c("gaussian","LOGNO"))) 
                    {
                    #subworkresp <- crossprod(bigZ, diag(x=new_weights))%*%c(new_workresp - X%*%new_beta - all_Zb - offset) ## REMOVE SINCE USING C++
                    #do.bvec <- try(chol2inv(chol(as.matrix(crossprod(sqrt(new_weights)*bigZ)) + kronecker(diag(x=length_uniids[k]),penmat2))) %*% subworkresp, silent=TRUE) ## REMOVE SINCE USING C++
                    do.bvec <- newcoeffn(sqrt(new_weights)*bigZ, kronecker(diag(x=length_uniids[k]),penmat2), c(new_workresp - X %*% new_beta - all_Zb - offset)*sqrt(new_weights)) 

                    if(!inherits(do.bvec,"try-error")) 
                         new_bvec <- do.bvec
                    if(inherits(do.bvec,"try-error")) {
                         subworkresp <- crossprod(bigZ, diag(x=new_weights))%*%c(new_workresp - X %*% new_beta - all_Zb - offset)
                         new_bvec <- ginv(as.matrix(crossprod(sqrt(new_weights)*bigZ)) + kronecker(diag(x=length_uniids[k]),penmat2)) %*% subworkresp 
                         }
                    }
               new_b[[k]] <- matrix(new_bvec, nrow = length_uniids[k], ncol = num_ran[k], byrow = TRUE)		
               sel_nonzerob <- which(sqrt(colSums(new_b[[k]]^2)) > length_uniids[k]*conv.eps)
               if(length(sel_nonzerob) > 0) 
                    new_b[[k]][,-sel_nonzerob] <- 0
               if(length(sel_nonzerob) == 0) 
                    new_b[[k]][,] <- 0
                
                
                ## Update random effects covariance matrix 
                ## Iterative estimator based on maximizing the Laplace approximation 
                ## If multiple Z are applied, then the Laplace approximated estimator is done conditonally on all other random effects
                new_ranef_covmat[[k]] <- matrix(0, nrow = num_ran[k], ncol = num_ran[k])
                if(length(sel_nonzerob) == 0) 
                    next;

                if(family$family[1] %in% c("gaussian","LOGNO")) {
                    for(i in 1:length_uniids[k]) { 
                         sel.i <- which(id[[k]] == i)
                         new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] <- new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] + chol2inv(chol(crossprod(as.matrix(Z[[k]][sel.i,sel_nonzerob,drop=FALSE]))/phi + do.inv.ranef_covmat[sel_nonzerob,sel_nonzerob])) + tcrossprod(new_b[[k]][i,sel_nonzerob])
                         }
                    new_ranef_covmat[[k]] <- new_ranef_covmat[[k]]/length_uniids[k]
                    }

                if(!(family$family[1] %in% c("gaussian","LOGNO"))) {
                    for(i in 1:length_uniids[k]) { 
                         sel.i <- which(id[[k]] == i)
                         new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] <- new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] + chol2inv(chol(crossprod(sqrt(new_weights[sel.i])*as.matrix(Z[[k]][sel.i,sel_nonzerob,drop=FALSE])) + do.inv.ranef_covmat[sel_nonzerob,sel_nonzerob])) + tcrossprod(new_b[[k]][i,sel_nonzerob])
                         }
                    new_ranef_covmat[[k]] <- new_ranef_covmat[[k]]/length_uniids[k]
                    }			

# 			if(ran.cov.est == "bb") { ## Based on (1/n_k) sum_{i=1}^{n_k} b_{ik}%*%t(b_{ik})
# 				new_ranef_covmat[[k]] <- cov.wt(new_b[[k]], center = FALSE, method = "ML")$cov
# 				}
# 				
# 			if(ran.cov.est == "wbb") { ## Based on (1/n_k) sum_{i=1}^{n_k} w_{ik} b_{ik}%*%t(b_{ik}), where w_{ik} is based on cluster size
# 				clus.size <- unlist(as.vector(table(id[[k]])))
# 				new_ranef_covmat[[k]] <- cov.wt(new_b[[k]], center = FALSE, method = "ML", wt = clus.size)$cov
# 				}
                        
                if(length(sel_nonzerob) == 0) 
                    ranef_fullshrunk[k] <- 1
                rm(bigZ, new_bvec)
                }

                    
          ## Solve scale parameters if required
          all_Zb <- numeric(length(y)) 
          for(k2 in 1:length(id)) {
               all_Zb <- all_Zb + as.numeric(all_bigZ[[k2]] %*% c(t(new_b[[k2]]))) 
               }
          if(family$family[1] == "gaussian") 
               new_phi <- mean((y - X %*% new_beta - all_Zb - offset)^2)
          if(family$family[1] == "LOGNO") 
               new_phi <- mean((log(y) - X %*% new_beta - all_Zb - offset)^2)

          if(family$family[1] == "negative.binomial") {
               new_phi <- try(1/theta.ml(y = y, mu = family$linkinv(X %*% new_beta + all_Zb + offset), limit = 100), silent = TRUE)
               if(inherits(new_phi,"try-error")) 
                    new_phi <- phi
               }
                
          if(family$family[1] == "Gamma") {
               prof.logl <- function(a, y, mu) {
                    out <- -y*a/mu + (a-1)*log(y) - lgamma(a) - a*log(mu) + a*log(a)
                    return(sum(out)) 
                    }
               update.shape <- try(optimize(f = prof.logl, interval = c(1e-3,1e3), y = y, mu = family$linkinv(X %*% new_beta+all_Zb+offset), maximum = TRUE), silent = TRUE)
               if(inherits(update.shape,"try-error"))    
                    new_shape <- shape 
               if(!inherits(update.shape,"try-error")) 
                    new_shape <- update.shape$max
               }
                
          if(family$family[1] == "ZIP") {
               prof.logl <- function(a, y, mu) {
                    out <- sum(dZIP(x=y, mu=mu, sigma=a, log=TRUE)) 
                    return(out)
                    }
               update.zeroprob <- try(suppressWarnings(optimize(f = prof.logl, interval = c(1e-4,0.9999), y = y, mu = family$mu.linkinv(X %*% new_beta+all_Zb+offset), maximum = TRUE)), silent = TRUE)
               if(inherits(update.zeroprob,"try-error")) 
                    new_zeroprob <- zeroprob
               if(!inherits(update.zeroprob,"try-error")) 
                    new_zeroprob <- update.zeroprob$max
               }
                    

                    
          diff <- mean((new_beta-beta)^2) 
          for(k in 1:length(id)) 
               diff <- diff + sum((new_ranef_covmat[[k]][upper.tri(new_ranef_covmat[[k]],diag=TRUE)]-ranef_covmat[[k]][upper.tri(ranef_covmat[[k]],diag=TRUE)])^2)
          #for(k in 1:length(id)) { diff <- diff + sum((new_b[[k]]-b[[k]])^2) }
          if(control$trace) { 
               cat("Iteration:", counter,"   Error:", round(diff,4), "\n")
                    
               new_beta2 <- new_beta
               names(new_beta2) <- colnames(X)
               message("Fixed effects:")
               print(c(new_beta2)); #print(c(newbeta.l2))
               message("Covariance Matrices:")
               print(new_ranef_covmat); 
               rm(new_beta2) 
               } 

                    
          beta <- new_beta
          b <- new_b; 
          ranef_covmat <- new_ranef_covmat; 
          shape <- new_shape
          zeroprob <- new_zeroprob
          phi <- new_phi
          counter <- counter + 1;
          true.counter <- true.counter + 1; 
          if(counter < 5) 
               true.counter <- 0

            
          if(any(unlist(lapply(ranef_covmat, function(x) any(sqrt(diag(x)) > 30)) == TRUE)) || sum(abs(new_beta)) < 1e-3) {
               cat("Convergence issues. Restarting...\n")
               restart_counter <- restart_counter + 1
               if(init_fit_run == 0) {
                    get_init <- start.fixed(y = y, X = X, Z = Z, id = id, family = family, offset = offset, trial.size = trial.size)
                    start <- build.start.fit.rpql(fit = get_init, id = id, num.ran = num_ran, cov.groups = cov.groups)
                    init_fit_run <- 1
                    }
                         
               new_beta <- beta <- start$fixef
               new_b <- b <- vector("list",length(id))
               for(k in 1:length(b)) { 
                    new_b[[k]] <- b[[k]] <- matrix(rnorm(length_uniids[k]*num_ran[k], mean=0, sd = 0.2+0.2*(family$family!="poisson")),length_uniids[k],num_ran[k]) 
                    }
               new_ranef_covmat <- ranef_covmat <- vector("list",length(id))
               for(k in 1:length(ranef_covmat)) 
                    new_ranef_covmat[[k]] <- ranef_covmat[[k]] <- cov(new_b[[k]])

               new_phi <- phi <- new_shape <- shape <- new_zeroprob <- zeroprob <- 1
               if(family$family[1] == "gaussian") 
                    new_phi <- phi <- mean((y - X %*% beta - offset)^2) 
               if(family$family[1] == "LOGNO") 
                    new_phi <- phi <- mean((log(y) - X %*% beta - offset)^2) 
               if(family$family[1] == "negative.binomial") 
                    new_phi <- phi <- 0.01 
               if(family$family[1] == "Gamma") 
                    new_shape <- shape <- 1
               if(family$family[1] == "ZIP") 
                    new_zeroprob <- zeroprob <- 0.1
                         
               diff <- 1e4
               cw_rpqllogLik <- 0
               counter <- 0 
               }		

          if(restart_counter > control$restarts) 
               { 
               stop("Sorry, but convergence could not be reached within the allocated number of control$restarts. Consider increasing lambda[2].\n")
               return()
               }

        }
        ## DONE!

		
     sel_nonzerob <- sel.zero.b <- vector("list",length(id))
     for(k in 1:length(id)) { 
          sel_nonzerob[[k]] <- which(sqrt(colSums(new_b[[k]]^2)) > length_uniids[k]*conv.eps) 
          if(length(sel_nonzerob[[k]]) > 0) 
               new_b[[k]][,-sel_nonzerob[[k]]] <- 0
          }
	
	
     ## The PQL likelihood
     new_etas <- X %*% new_beta + offset
     for(k2 in 1:length(id)) 
          new_etas <- new_etas + as.numeric(all_bigZ[[k2]] %*% c(t(new_b[[k2]])))

     if(family$family[1] == "gaussian") 
          pqllogLik <- sum(dnorm(y, mean=new_etas, sd=sqrt(new_phi), log=TRUE)) 
     if(family$family[1] == "Gamma") 
          pqllogLik <- sum(dgamma(y, shape=new_shape, scale=family$linkinv(new_etas)/new_shape, log=TRUE)) 
     if(family$family[1] == "binomial") 
          pqllogLik <- sum(dbinom(y, size=trial.size, prob=family$linkinv(new_etas), log=TRUE)) 
     if(family$family[1] == "poisson") 
          pqllogLik <- sum(dpois(y, lambda=family$linkinv(new_etas), log=TRUE)) 
     if(family$family[1] == "negative.binomial") 
          pqllogLik <- sum(dnbinom(y, mu=family$linkinv(new_etas), size=1/new_phi, log=TRUE))
     if(family$family[1] == "LOGNO") 
          pqllogLik <- sum(dlnorm(y, meanlog=new_etas, sdlog=sqrt(new_phi), log=TRUE))
     if(family$family[1] == "ZIP") 
          pqllogLik <- sum(dZIP(y, mu=family$mu.linkinv(new_etas), sigma=new_zeroprob, log=TRUE))
     loglik1 <- pqllogLik
    
     for(k in 1:length(id)) { 
          if(length(sel_nonzerob[[k]]) > 0) 
               pqllogLik <- pqllogLik - 0.5*sum(rowSums(new_b[[k]][,sel_nonzerob[[k]]]*(new_b[[k]][,sel_nonzerob[[k]]]%*%solve(new_ranef_covmat[[k]][sel_nonzerob[[k]],sel_nonzerob[[k]]])))) 
          }
     names(sel_nonzerob) <- names(id)
               
               
     ## Information Criterion...argh!
     aic.v1 <- -2*loglik1 + 2*(sum(new_beta != 0) + sum(sapply(new_b, function(x) sum(x!=0))))
     bic.v2 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + sum(log(length_uniids)*sapply(new_ranef_covmat, function(x) sum(x[upper.tri(x,diag=T)]!=0))) ## Number of random effects to penalize is just number of non-zero elements in covariance matrix
     bic.v3 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + sum(log(length(y))*sapply(new_ranef_covmat, function(x) sum(x[upper.tri(x,diag=T)]!=0))) ## Number of random effects to penalize is just number of non-zero elements in covariance matrix
     hic.v1 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + 2*sum(sapply(new_b, function(x) sum(x!=0)))
     hic.v2 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + sum(sapply(new_b, function(x) sum(x!=0)))
     hic.v3 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + 0.5*sum(sapply(new_b, function(x) sum(x!=0)))
	
	
     dof <- sum(new_beta != 0) + sum(sapply(new_b, function(x) sum(x!=0)))
     ics <- c(dof, aic.v1, bic.v2, bic.v3, hic.v1, hic.v2, hic.v3) 
     names(ics) = c(
          "# of estimated parameters",
          "AIC: 2*sum(beta!=0) + 2*sum(b!=0)",
          "BIC: log(nrow(X))*sum(beta!=0) + log(number of clusters)*sum(ranef_covmat!=0)",
          "BIC: log(nrow(X))*sum(beta!=0) + log(nrow(X))*sum(ranef_covmat!=0)",
          "Hybrid IC: log(nrow(X))*sum(beta!=0) + 2*sum(b!=0)",
          "Hybrid IC: log(nrow(X))*sum(beta!=0) + sum(b!=0)",
          "Hybrid IC: log(nrow(X))*sum(beta!=0) + 0.5*sum(b!=0)")

               
     ## Garnishing
     out_list <- list(fixef = as.vector(new_beta), ranef = new_b, ran.cov = new_ranef_covmat, pql.logLik = pqllogLik, logLik1 = loglik1, 
          phi = new_phi, shape = new_shape, zeroprob = new_zeroprob, family = family, n = length(y), 
          trial.size = trial.size, id = id, lambda = lambda, pen.weights = pen.weights, pen.type = pen.type, control = control,
          ics = ics, nonzero.fixef = which(as.vector(new_beta)!=0), nonzero.ranef = sel_nonzerob)

     if(save.data) { 
          out_list$y <- y
          out_list$X <- X
          outlist$Z <- Z
          out_list$offset <- offset
          message("Please note the column names in Z have been overwritten for easier reference in the estimation algorithm. Apologies in advance") 
          }

		
     if(hybrid.est) {
          if(family$family[1] %in% c("ZIP","LOGNO")) { 
               message("Hybrid estimation not possible as lme4 does not do this. Sorry! Maybe try glmmTMB?")
               out_list$hybrid <- NULL 
               }
                    
          else {	
               if(control$trace) 
                    message("Performing hybrid estimation with lme4...might take a bit of time, and the labels for the parameters will be wrong. Apologies in advance!")

               make_dat <- data.frame(y, do.call(cbind,Z), X, do.call(cbind,id))
               restring1 <- NULL
               if(length(out_list$nonzero.fixef) > 0) 
                    restring1 <- paste(paste(colnames(X)[out_list$nonzero.fixef],collapse="+"),"-1")
               for(k in 1:length(id)) { 
                    if(length(out_list$nonzero.ranef[[k]])) restring1 <- paste(restring1, "+ (",paste0(colnames(Z[[k]])[out_list$nonzero.ranef[[k]]],collapse="+"),"-1|",names(id)[k],")")
                    }
               restring1 <- reformulate(restring1, response = "y")
    # 		print(restring1)
               if(family$family[1] == "gaussian") 
                    final_fit <- suppressWarnings(lmer(restring1, REML = TRUE, offset = offset, data = make_dat))
               if(family$family[1] != "gaussian" & family$family[1] != "negative.binomial") 
                    final_fit <- suppressWarnings(glmer(restring1, family = family, offset = offset, data = make_dat))
               if(family$family[1] == "negative.binomial") 
                    final_fit <- suppressWarnings(glmer.nb(restring1, family = family, offset = offset, data = make_dat))

               out_list$hybrid <- final_fit ## Access covariance matrix as VarCorr(final_fit)[[1]][1:nrow(VarCorr(final_fit)[[1]]),1:nrow(VarCorr(final_fit)[[1]])]
                    
               rm(make_dat)
               }
          }	
    

     names(out_list$fixef) <- colnames(X)
     names(out_list$ranef) <- names(out_list$ran.cov) <- names(id)
     for(k in 1:length(id)) { 
          rownames(out_list$ranef[[k]]) <- 1:length_uniids[k]; colnames(out_list$ranef[[k]]) <- colnames(Z[[k]])
          rownames(out_list$ran.cov[[k]]) <- colnames(out_list$ran.cov[[k]]) <- colnames(Z[[k]]) 
          }
            
     class(out_list) <- "rpql"
     out_list$call <- match.call()
    
     return(out_list)
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
	
	
