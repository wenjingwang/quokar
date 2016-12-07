#'Outlier Dignostic for Quantile Regression Based on Bayesian Estimation
#'@param y dependent variable in quantile regression
#'
#'@param x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix
#'
#'@param tau quantile
#'
#'@param M the iteration frequancy for MCMC used in Baysian Estimation
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
#'The method cook.distance, qfunction.distance
#'
#'
#'@author Wenjing Wang \email{wenjingwang1990@gmail.com}
#'@keywords quantile regression outlier diagnostics
#'@references Benites L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#' Diagnostics for Quantile Regression Using the Asymmetric Laplace
#' Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#'Hawkins D M, Bradu D, Kass G V.(1984)``Location of several outliers in
#'multiple-regression data using elemental sets. \emph{Technometrics},
#'26(3), 197-208.
#'
#'Koenker R, Bassett Jr G.(1978)`` Regression quantiles,
#'\emph{Econometrica}, 1, 33-50.
#'
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'Kozumi H, Kobayashi G.(2011)``Gibbs sampling methods for Bayesian
#'quantile regression,\emph{Journal of statistical computation and
#'simulation}, 81(11), 1565-1578.
#'
#'
#'@export
#'
#'



qrod_bayes <- function(y, x, tau, M, method = c("bayes.prob","bayes.kl")){
  method <- match.arg(method)
  if(method == "bayes.prob"){
    result <- bayesProb(y, x, tau, M)
  }else if(method == "bayes.kl"){
    result <- bayesKL(y, x, tau, M)
  }else if(method %in% c("bayes.prob","bayes.kl") == FALSE)
    warning("Method should be one of 'bayes.prob', 'bayes.kl'")
  return(result)
}

