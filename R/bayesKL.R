#'Kullback-Leibler divergence for each observation in Baysian quantile regression model
#'@param y dependent variable in quantile regression
#'
#'@param x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix
#'
#'@param tau quantile
#'
#'@param M the iteration frequancy for MCMC used in Baysian Estimation
#'
#'@param burn burned MCMC draw
#'
#'@details
#'Method to address the differences between the posterior distributions
#'from the distinct latent variables in the model, we suggest the use of the Kullback-
#'Leibler divergence as a more precise method of measuring the distance between those
#'latent variables in the Bayesian quantile
#'regression framework. In this posterior information, the divergence is defined as
#'
#'\deqn{K(f_{i}, f_{j}) = \int log(\frac{f_{i}(x)}{f_{j}{(x)}})f_{i}(x)dx}
#'
#'where \eqn{f_{i}} could be the posterior conditional distribution of \eqn{v_{i}}
#'and \eqn{f_{j}} the poserior conditional distribution of \eqn{v_{j}}. We should average
#'this divergence for one observation based on the distance from all others, i.e,
#'
#'\deqn{KL(f_{i})=\frac{1}{n-1}\sum{K(f_{i}, f_{j})}}
#'
#'We expect that when an observation presents a higher value for this divergence,
#'it should also present a high probability value of being an outlier. Based on
#'the MCMC draws from the posterior of each latent vaiable, we estimate the densities
#'using a normal kernel and we compute the integral using the trapezoidal rule.
#'
#'More details please refer to the paper in references
#'@references
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'@seealso \code{bayesProb}
#'
#'
#'
#'
#'

trapz <- function(x,y){
  idx = 2:length(x)
  return (as.double((x[idx]-x[idx-1]) %*% (y[idx]+y[idx-1]))/2)
}
bayesKL <- function(y, x, tau, M, burn){
  x <- cbind(1, x)
  n <- length(y)
  t <- M - burn
  coefs <- bayesQR(y ~ x, quantile = tau, alasso = FALSE,
                   ndraw = M, prior = NULL)
  beta <- coefs[[1]]$betadraw[(burn+1):M, ]
  sigma <- coefs[[1]]$sigmadraw[(burn+1):M]
  taup2 <- (2/(tau * (1 - tau)))
  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  v <- matrix(0, nrow = n, ncol = t)
  for (i in 1:n) {
    for (j in 1:t) {
      param1 <- 1/2
      param2 <- (y[i] - x[i, ] %*% t(beta[j,]))^2/(taup2*sigma[j])
      param3 <- 2/sigma[j] + theta^2/(taup2*sigma[j])
      v[i, j] <- rgig(1, param1, param2, param3)
    }
  }
  KLS <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      hi <- stats::density(v[i,], kernel = "gaussian")$bw
      hj <- stats::density(v[j,], kernel = "gaussian")$bw
      func1 <- function(x, y, h){
        mid <- ((x-y)/h)^2
        fun1 <- sum((1/sqrt(2*pi))*exp(-mid/2))/(t*h)
        return(fun1)
      }
      func2 <- function(x, ...){
        fun2 <- log(func1(x, v[i,], hi)/func1(x, v[j,], hj))*func1(x, v[i,], hi)
        return(fun2)
      }
      upper_x <- max(max(v[i,],v[j,]))
      lower_x <- min(min(v[i,],v[j,]))
      xx <- seq(lower_x,upper_x,length.out = 100)
      yy <- rep(0,length(xx))
      for(k in 1:length(xx)){
        yy[k] <- func2(xx[k])
      }
      KLS[i,j] <- trapz(xx,yy)
    }
  }
  KLD <- rowSums(KLS)/(n-1)
  return(KLD)
}

