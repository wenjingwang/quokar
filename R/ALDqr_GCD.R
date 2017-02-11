#'Calculating generalized cook distance of the MLE estimation of quantile regression using asymmetric laplace distribution
#'@param y Dependent variable in quantile regression. Note that: we suppose
#'y follows asymmetric laplace distribution.
#'
#'@param x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix which including
#'the intercept. That means, if the dimension of independent
#'variables is p and the sample size is n, x is a n times p+1
#'matrix with the first column is one.
#'
#'@param tau quantile
#'
#'@param error The EM algorithm accuracy of error used in
#' MLE estimation
#'
#'@param iter the iteration frequancy for EM algorithm used
#' in MLE estimation
#'
#'
#'@details
#'Case-deletion is a classical approach to study the effects
#' of dropping the
#'\eqn{i}th case from the data set. Thus, the complete-data log-likelihhod
#'function based on the data with \eqn{i}th cse deleted with be denoted by
#'\eqn{l_{c}(\theta|y_{c(i)})}. Let \eqn{\hat{\theta_{p(i)}} = (\hat{\beta^{'}_{p(i)}},
#'\hat{\sigma^{2}}_{(i)})^{'}} be the maximizer of the function
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
#'
#'@references
#'Benites L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#'Diagnostics for Quantile Regression Using the Asymmetric Laplace
#'Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#'@seealso \code{ALDqr_QD}
#'
#'
ALDqr_GCD <- function(y, x, tau, error, iter)
{
  n <- length(y)
  p <- ncol(x)
  theta <- ALDqr::EM.qr(y, x, tau, error, iter)$theta
  beta_qr <- theta[1:p, ]
  sigma_qr <- theta[p+1]
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
  gamma2 <- (2 + thep^2/taup2)/sigma_qr
  muc <- y - x %*% beta_qr
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/
      (besselK(sqrt(delta2* gamma2), 0.5))*(sqrt(delta2 / gamma2))^(-1)

  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2              *gamma2), 0.5)) * (sqrt(delta2 / gamma2))
  E1 <- matrix(0, nrow = p, ncol = n)
  for(i in 1:n){
      suma2 <- x[-i,] * c(vchpN[-i] * (y[-i] - x[-i,] %*% beta_qr) - thep)
      E1[,i] <- apply(suma2, 2, sum)/(taup2)
  }
  E2 <- 1: n %>%
      map(function(i) {
      muc_i <- y[-i] - x[-i, ]%*%beta_qr
      sum(3*sigma_qr - (vchpN[-i] * muc_i^2 -
                          2 * muc_i * thep + vchp1[-i] *(thep^2 + 2 * taup2))/taup2)
    })
  E2 <- simplify2array(E2)
  Q1_beta <- E1/sigma_qr
  Q1_sigma <- -E2/(2*sigma_qr^2)
  xM <- c(sqrt(vchpN)) * x
  suma1 <- t(xM) %*% (xM)
  Q2_beta <- -(suma1)/(sigma_qr * taup2)
  Q2_sigma <- 3/(2*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +
                                        vchp1*(thep^2 + 2*taup2)))/(sigma_qr^3*taup2)
  GCD_beta <- 1:n %>%
    map(function(i){
      c(Q1_beta[,i]) %*% solve(-Q2_beta) %*% matrix(Q1_beta[,i], ncol = 1)
    })
  GCD_beta <- simplify2array(GCD_beta)
  GCD_sigma <- 1:n %>%
    map(function(i) {
      Q1_sigma[i]*solve(-Q2_sigma)*Q1_sigma[i]
    })
  GCD_sigma <- simplify2array(GCD_sigma)
  GCD <- GCD_beta + GCD_sigma
  return(GCD)
}
