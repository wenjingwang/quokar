#'Calculating Kullback-Leibler divergence vased on the bayes estimating procedure of
#'of quantile regression with asymmetric laplace distribution
#'@param y dependent variable in quantile regression
#'
#'@param x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix
#'
#'@param tau quantile
#'
#'@param M the iteration frequancy for MCMC used in Baysian Estimation
#'
#'@details
#'Method to address the differences between the posterior distributions
#'from the distinct latent variables in the model, we suggest the use of the Kullback-
#'Leibler divergence proposed by Kullback and Leibler(1951), as a more precise method
#'of measuring the distance between those latent variables in the Bayesian quantile
#'regression framework. In this posterior information, the divergence is defined as
#'
#'\deqn{K(f_{i}, f_{j}) = \int log(\frac{f_{i}(x)}{f_{j}{(x)}})f_{i}(x)dx}
#'
#'where \eqn{f_{i}} could be the posterior conditional distribution of \eqn{v_{i}}
#'and \eqn{f_{j}} the poserior conditional distribution of \eqn{v_{j}}. Similar to
#'the probability proposal in the previous subsection, we should average this
#'divergence for one observation based on the distance from all others, i.e,
#'
#'\deqn{KL(f_{i})=\frac{1}{n-1}\sum{K(f_{i}, f_{j})}}
#'
#'We expect that when an observation presents a higher value for this divergence,
#'it should also present a high probability value of being an outlier. Based on
#'the MCMC draws from the posterior of each latent vaiable, we estimate the densities
#'using a normal kernel and we compute the integral using the trapezoidal rule.
#'
#'@references
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'@seealso \code{bayesProb}
#'
bayesKL <- function(y, x, tau, M){
  coefs <- bayesQR::bayesQR(y ~ x, quantile = tau, alasso = FALSE, ndraw = M)
  beta <- summary(coefs)[[1]]$betadraw[, 1]
  sigma <- summary(coefs)[[1]]$sigmadraw[, 1]
  taup2 <- (2/(tau * (1 - tau)))
  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  n <- length(y)
  x <- cbind(1, x)
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
