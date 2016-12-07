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
bayesKL <- function(y, x, tau, M){
  coefs <- bayesQR::bayesQR(y ~ x, quantile = tau, alasso = FALSE, ndraw = M)
  beta <- summary(coefs)[[1]]$betadraw[, 1]
  sigma <- summary(coefs)[[1]]$sigmadraw[, 1]
  taup2 <- (2/(tau * (1 - tau)))
  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  n <- length(y)

  delta <- sqrt((y - x %*% beta)^2/(sigma * taup2))
  gamma <- sqrt(theta^2/(sigma*taup2) + 2/sigma)

  nu_dens <- matrix(0, nrow = M, ncol = n)

  for(i in 1:n) {
    nu_dens[,i] <- HyperbolicDist::rgig(M, 1/2, delta[i], gamma)
  }

  KLS <- matrix(0, nrow = n, ncol = n-1)

  for(i in 1:n) {
    nu1 <- nu_dens[, -i]
    hi <- stats::density(nu_dens[, i], kernel = "gaussian")$bw
    upper_x <- max(nu_dens[, i])
    lower_x <- min(nu_dens[, i])
    for(j in 1:(n-1)) {
      hj <- stats::density(nu1[, j], kernel = "gaussian")$bw
      func <- function(xx) {
        log((1/(M*sqrt(2*pi))*sum(exp(-1/2*((xx-nu_dens[,i])/hi)^2)))/
              (1/(M*sqrt(2*pi))*sum(exp(-1/2*((xx-nu1[,j])/hj)^2))))*
          (1/(M*sqrt(2*pi))*sum(exp(-1/2*((xx-nu_dens[,i])/hi)^2)))
      }
      KLS[i,j] <- (upper_x - lower_x)*(func(upper_x)+func(lower_x))/2
    }
  }

  KLD <- apply(KLS, 1, mean)
  for(i in 1:n) {
    if (KLD[i] == Inf) {KLD[i] = 1}
  }
  return(KLD)
}
