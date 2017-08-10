#'@title Outlier Dignostic for Quantile Regression with Asymmetric
#'Laplace Distribution
#'@param y Dependent variable in quantile regression
#'
#'@param x Indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix which including
#'the intercept. That means, if the dimension of independent
#'variables is p and the sample size is n, x is a n times p+1
#'matrix with the first column being one.
#'
#'@param tau quantile
#'
#'@param error The EM algorithm accuracy of error used in MLE
#'estimation
#'
#'@param iter the iteration frequancy for EM algorithm used in
#'MLE estimation
#'
#'@param method the diagnostic method for outlier detection
#'
#'@return Generalized Cook's distance or Q-function distance for
#' multiple quantiles
#'
#'@description
#'This function cacluate the generalized cook distance and q
#'function distance of quantile regression model with
#'asymmetric laplace distribution.
#'
#'@details please refer to the reference paper
#'
#'
#'
#'
#'
qrod_mle <- function(y, x, tau, error, iter,
                 method = c("cook.distance", "qfunction")){
  method <- match.arg(method)
  ntau <- length(tau)
  n <- length(y)
  distance <- list()
  if(!(method %in% c("cook.distance", "qfunction"))){
    stop("Method should be 'cook.distance' or 'qfunction'")
  }else if(method == "cook.distance"){
      for(i in 1:ntau){
          distance[[i]] <- ALDqr_GCD(y, x, tau[i], error, iter)
      }
    }else if(method == "qfunction"){
        for(i in 1:ntau){
            distance[[i]] <- ALDqr_QD(y, x, tau[i], error, iter)
      }
    }
  return(distance)
}





