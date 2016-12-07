#'Outlier Dignostic for quantile regression
#'@para y dependent variable in quantile regression
#'
#'@para x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix which including
#'the intercept. That means, if the dimension of independent
#'variables is p and the sample size is n, x is a (n \times (p+1))
#'matrix with the first column is 1.
#'
#'@para error The EM algorithm accuracy of error used in MLE estimation
#'
#'@para iter the iteration frequancy for EM algorithm used in MLE estimation
#'
#'@M the iteration frequancy for MCMC used in Baysian Estimation
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
#'Hawkins D M, Bradu D, Kass G V.(1984)``Location of several outliers in
#'multiple-regression data using elemental sets. \emph{Technometrics},
#'\bold{26(3)}, 197-208.
#'
#'Koenker R, Bassett Jr G.(1978)`` Regression quantiles,
#'\emph{Econometrica},\bold{1}, 33-50.
#'
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'Kozumi H, Kobayashi G.(2011)``Gibbs sampling methods for Bayesian
#'quantile regression,\emph{Journal of statistical computation and
#'simulation}, \bold{81(11)}, 1565-1578.
#'
#'@examples
#'data(ais, package = `ALDqr`)
#'
#'
#'
#'

qrod <- function(y, x, tau, error = NULL, iter = NULL,
                 method = c("cook.distance", "qfunciton", "bayes.probability",
                            "bayes.KLD"), M = NULL){
  method <- match.arg(method)
    if(method == "cook.distance"){
      distance <- ALDqr_GCD(y, x, tau, error, iter)
    }else if(method == "qfunction"){
      distance <- ALDqr_QD(y, x = NULL, tau = NULL, error, inter)
    }else if(method == "bayes.probability"){
      distance <- bayesProb(y, x, beta, sigma, tau, M)
    }else if(method == "bayes.KLD"){
      distance <- bayesKLD(y, x, beta, sigma, tau, M)
    }else if(method %in% c("cook.distance", "qfunction", "bayes.probability",
                            "bayes.KLD") == FALSE)
      warning("Method should be one of 'cook.distance', 'qfunction',
              'bayes.probability', 'bayes.KLD'")
  returen(distance)
}

