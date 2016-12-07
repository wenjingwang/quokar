#'Outlier Dignostic for Quantile Regression Based on MLE Estimation
#'@param y dependent variable in quantile regression
#'
#'@param x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix which including
#'the intercept. That means, if the dimension of independent
#'variables is p and the sample size is n, x is a n times p+1
#'matrix with the first column is one.
#'
#'@param tau quantile
#'
#'@param error The EM algorithm accuracy of error used in MLE estimation
#'
#'@param iter the iteration frequancy for EM algorithm used in MLE estimation
#'
#'@param method the diagnostic method for outlier detection
#'
#'@description
#'This group of function is used to compute diagnositcs for a
#'quantile regression model based on the building blocks returned
#'by.
#'
#'@details
#'The primary function is which returns either
#'a list or data frame of influence mesures depending on whether
#'or if only one aspect of the model is selected
#'then a list with Generalized Cook's distance,Q-function Distance,
#'Baysian Method or ES method.
#'
#'This function is used to detect outlier in quantile regression
#'level. It use code to fit the models for each quantile.
#'
#'The method cook.distance, qfunction.distance,
#'baysian.probability, baysian.distribution and ES
#'can be used for direct computation of the corresponding diagnostic
#'quantities from an object of class case_delete
#'
#'@author Wenjing Wang \email{wenjingwang1990@gmail.com}
#'@keywords quantile regression outlier diagnostics
#'@references Benites L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#' Diagnostics for Quantile Regression Using the Asymmetric Laplace
#' Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#'
#'Koenker R, Bassett Jr G.(1978)`` Regression quantiles,
#'\emph{Econometrica},\bold{1}, 33-50.
#'
#'
#'@export
#'
#'

qrod_mle <- function(y, x, tau, error, iter,
                 method = c("cook.distance", "qfunciton")){
  method <- match.arg(method)
    if(method == "cook.distance"){
      distance <- ALDqr_GCD(y, x, tau, error, iter)
    }else if(method == "qfunction"){
      distance <- ALDqr_QD(y, x, tau, error, iter)
    }else if(method %in% c("cook.distance", "qfunction") == FALSE)
      warning("Method should be one of 'cook.distance', 'qfunction'")
  return(distance)
}

