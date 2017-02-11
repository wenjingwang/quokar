#'Calculating the case-deletion coefficience of quantile regression
#'@param y Response variable in quantile regression model
#'
#'@param x Predictors in quantile regression model.
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
#'@importFrom purrr %>% map
#'
#'
#'
ALDqr_case_deletion <- function(y, x, tau, error, iter)
{
  n <- length(y)
  p <- ncol(x)
  qr <- ALDqr::EM.qr(y,x,tau,error,iter)
  beta_qr <- qr$theta[1:p,]
  sigma_qr <- qr$theta[p+1]
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
  gamma2 <- (2 + thep^2/taup2)/sigma_qr
  muc <- y - x %*% beta_qr
  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                      gamma2), 0.5)) * (sqrt(delta2 / gamma2))
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                      gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)
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
  beta_i <- matrix(0, nrow=p, ncol = n)
  for(i in 1:n){
    beta_i[,i] <- beta_qr + taup2*solve(suma1)%*% E1[,i]
  }
  sigma_i2 <- 1:n %>%
    map(function(i) sigma_qr^2 - solve(Q2_sigma)*E2[i]/(2*sigma_qr^2))
  sigma_i <- sqrt(simplify2array(sigma_i2))
  theta_i <- list(beta_i = beta_i, sigma_i = sigma_i)
  return(theta_i)
}

