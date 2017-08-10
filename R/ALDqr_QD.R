#'Q-function distance for each observation in quantile regression model
#'@param y Dependent variable in quantile regression. Note that: we suppose
#'y follows asymmetric laplace distribution.
#'
#'@param x Indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix which including
#'the intercept. That means, if the dimension of independent
#'variables is p and the sample size is n, x is a n times p+1
#'matrix with the first column is one.
#'
#'@param tau Quantile
#'
#'@param error The EM algorithm accuracy of error used in MLE estimation
#'
#'@param iter The iteration frequancy for EM algorithm used in MLE estimation
#'
#'@details
#'Measure of the influence of the \eqn{i}th case is the following Q-distance
#'function, similar to the likelihood distance \eqn{LD_{i}} (Cook and Weisberg, 1982),
#'defined as
#'
#'\deqn{QD_{i} = 2{Q(\hat{\theta}|\hat{\theta})-Q(\hat{\theta_{(i)}})}}
#'
#'@references
#'Benites L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#'Diagnostics for Quantile Regression Using the Asymmetric Laplace
#'Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#'@author Wenjing Wang \email{wenjingwangr@gmail.com}
#'@seealso \code{ALDqr_GCD}

ALDqr_QD <- function(y, x, tau, error, iter)
{
  p <- ncol(x)
  n <- length(y)
  theta_all <- ALDqr::EM.qr(y, x, tau, error, iter)$theta
  beta_all <- theta_all[1:p, ]
  sigma_all <- theta_all[p+1]
  beta_i <- ALDqr_case_deletion(y, x, tau, error, iter)$beta_i
  sigma_i <-  ALDqr_case_deletion(y, x, tau, error, iter)$sigma_i
  Q_function <- function(beta_qr, sigma_qr, tau){
    n <- length(y)
    taup2 <- (2/(tau * (1 - tau)))
    thep <- (1 - 2 * tau) / (tau * (1 - tau))
    delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
    gamma2 <- (2 + thep^2/taup2)/sigma_qr
    muc <- y - x %*% beta_qr
    vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                    gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)
    vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                    gamma2), 0.5)) * (sqrt(delta2 / gamma2))
    Q <- (-3*log(sigma_qr)/2)*n - sum(vchpN * muc^2 - 2 * muc * thep +
                vchp1 *(thep^2 + 2 * taup2))/(2 * sigma_qr * taup2)
    return(Q)
  }
  Q_all <- Q_function(beta_all, sigma_all, tau)
  Q_i <- rep(0, n)
  for(i in 1:n){
    Q_i[i] <- Q_function(beta_i[,i], sigma_i[i], tau)
  }
  QD <- rep(0, n)
  for(i in 1:n){
    QD[i] <- 2*(Q_all - Q_i[i])
  }
  return(QD)
}
