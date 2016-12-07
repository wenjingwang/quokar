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

bayesProb <- function(y, x, tau, M){
   coefs <- bayesQR::bayesQR(y ~ x, quantile = tau, alasso = FALSE,
                    ndraw = M, prior = NULL)
   beta <- summary(coefs)[[1]]$betadraw[, 1]
   sigma <- summary(coefs)[[1]]$sigmadraw[, 1]
   taup2 <- (2/(tau * (1 - tau)))
   theta <- (1 - 2 * tau) / (tau * (1 - tau))
   delta <- sqrt((y - x %*% beta)^2/(sigma * taup2))
   gamma <- sqrt(theta^2/(sigma*taup2) + 2/sigma)
   n = length(y)
   nu_dens <- matrix(0, nrow = M, ncol = n)
   for(i in 1:n) {
       nu_dens[,i] <- HyperbolicDist::rgig(M, 1/2, delta[i], gamma)
   }
   A <- matrix(0, nrow = n, ncol = n-1)
   for(i in 1:n) {
       probs <- 1:M/M
       nu1 <- nu_dens[, -i]
       for(j in 1:(n-1)) {
       A[i,j] <- 1/M*sum(stats::quantile(nu_dens[,i], probs = probs) > max(nu1[, j]))
       }
   }
prob <- apply(A, 1, mean)
return(prob)
}
