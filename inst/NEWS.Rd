\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN rpql VERSION 0.8}{
This list below is written in terms of what is deemed to be most important to least important changes =P 

  \itemize{
     \item A new argument \code{lasso.lambda.scale} has been included as part of \code{control} in the main \code{rpql} function. This scales the penalty on the group random effects, when the lasso or adaptive lasso penalty is used, by \code{1/sqrt(no. of clusters)} so that it is on the same scale as the fixed effects penalty; Thanks to Sarbesh Pandeya for this suggestion. Defaults to TRUE.
     
     \item Modified some of the examples in response to the above change.
     
     \item Fixed a bug issue when \code{trial.size > 1}.
     }
}
