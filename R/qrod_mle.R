#'@title Outlier Dignostic for Quantile Regression with Asymmetric Laplace Distribution
#'@param y Dependent variable in quantile regression
#'
#'@param x Indepdent variables in quantile regression.
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
#'
#'This function cacluate the generalized cook distance and q function distance
#'of quantile regression model with asymmetric laplace distribution.
#'
#'@details
#'
#'Case-deletion is a classical approach to study the effects of dropping the
#'\eqn{i}th case from the data set. Thus, the complete-data log-likelihhod
#'function based on the data with \eqn{i}th cse deleted with be denoted by
#'\eqn{l_{c}(\theta|y_{c(i)})}. Let \eqn{\hat{\theta_{p(i)}} = (\hat{\beta^{'}_{p(i)}},
#'\hat{\sigam^{2}}_{(i)})^{'}} be the maximizer of the function
#'
#'\deqn{Q_{(i)}(\theta|\hat{\theta})=E_{\hat{\theta}}[l_{c}(\theta|Y_{c(i)})|y]}
#'
#'To assess the influence of the \eqn{i}th case on the EM estimate \eqn{\hat{\theta}},
#'we compare \eqn{\hat{\theta_(i)}} and \eqn{\hat{\theta}}, and if \eqn{\hat{\theta_(i)}}
#'is far from \eqn{\hat{\theta_(i)}} in some sense, then the \eqn{i}th case is regarded
#'as influential. Based on the metric for measuring the distance between \eqn{\hat{\theta_(i)}}
#'and \eqn{\hat{\theta}} proposed by Zhu et al.(2001), we consider here the following
#'generalized Cook distance:
#'
#'\deqn{GD_{i} = (\hat{\theta_{(i)}}-\hat{\theta{i}})^{'}{-Q(\hat{\theta}|\hat{\theta})}
#'(\hat{\theta_{(i)}}-\hat{\theta{i}})}
#'
#'Another measure of the influence of the \eqn{i}th case is the following Q-distance
#'function, similar to the likelihood distance \eqn{LD_{i}} (Cook and Weisberg, 1982),
#'defined as
#'
#'\deqn{QD_{i} = 2{Q(\hat{\theta}|\hat{\theta})-Q(\hat{\theta_{(i)}})}}
#'
#'@author Wenjing Wang \email{wenjingwang1990@gmail.com}
#'@keywords quantile regression outlier diagnostics
#'
#'@references
#'Sánchez B L, Lachos H V, Labra V F.(2013)``Likelihood based inference for quantile
#'regression using the asymmetric Laplace distribution,\emph{
#'Journal of Statistical Computation and Simulation}, 81: 1565-1578.
#'
#'L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#'Diagnostics for Quantile Regression Using the Asymmetric Laplace
#'Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#'
#'Koenker R, Bassett Jr G.(1978)`` Regression quantiles,
#'\emph{Econometrica},\bold{1}, 33-50.
#'
#'@seealso \code{covMCD} from package `robustbase`
#'@examples
#'data(ais)
#'y <- ais$BMI
#'sexInd <- (ais$Sex == 1) + 0
#'x <- cbind(1, ais$LBM, ais$sexInd)
#'qrod_mle(y, x, tau = 0.1, error = 1e-10, iter = 2000, method = "cook.distance")
#'qrod_mle(y, x, tau = 0.1, error = 1e-10, iter = 2000, method = "qfunction")
#'
#'
#'
#'@export
#'
#'
qrod_mle <- function(y, x, tau, error, iter,
                 method = c("cook.distance", "qfunciton")){
  method <- match.arg(method)
  if(!(method %in% c("cook.distance", "qfunction"))){
    stop("Method should be 'cook.distance' or 'qfunction'")
    }else if(method == "cook.distance"){
      distance <- ALDqr_GCD(y, x, tau, error, iter)
      case_distance <- data.frame(case = 1:length(y), distance = distance)
    }else if(method == "qfunction"){
      distance <- ALDqr_QD(y, x, tau, error, iter)
      case_distance <- data.frame(case = 1:length(y), distance = distance)
    }
  return(case_distance)
}

