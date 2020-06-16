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
     for(k in 1:length(x$ranef)) { 
          message(names(x$id)[k], ":")
          print(diag(x$ran.cov[[k]]))
          message() 
          } 

    if(x$family$family[1] %in% c("gaussian","lognormal")) 
        message("Variance par ameter: ", x$phi, "\n")
    if(x$family$family[1] %in% c("negative.binomial")) 
        message("overdispersion parameter (V = mu + phi*mu^2): ", x$phi, "\n")
    if((x$family$family[1] == "Gamma")) 
        message("Shape parameter (V = mu^2/shape): ", x$shape, "\n")
    if((x$family$family[1] == "ZIP")) 
        message("Probability of structural zero: ", x$zeroprob, "\n") 
    }	
	
	
summary.rpql <- function(object, ...) {
     gather_output <- list(call = object$call, fixef = round(object$fixef,3), ranef = lapply(object$ranef,round,3), 
          ran.cov = lapply(object$ran.cov,round,3), logLik = round(object$logLik,3), family = object$family)
    
    if(object$family$family[1] %in% c("gaussian","lognormal","negative.binomial")) 
        gather_output$phi <- round(object$phi,3)
    if((object$family$family[1] == "Gamma")) 
        gather_output$shape <- round(object$shape,3)
    if((object$family$family[1] == "ZIP")) 
        gather_output$zeroprob <- round(object$zeroprob,3)

    gather_output$pen.type <- object$pen.type
    gather_output$lambda <- object$lambda
    gather_output$ics <- object$ics
    gather_output$id <- object$id 
    gather_output$nonzero.fixef <- object$nonzero.fixef
    gather_output$nonzero.ranef <- object$nonzero.ranef
    gather_output$ics <- object$ics 
    
    class(gather_output) <- "summary.rpql"
    gather_output 
    }	
	
	
	
