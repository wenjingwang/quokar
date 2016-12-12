#'Calculating mean posterior probability based on the bayes estimating procedure of
#'quantile regression with asymmetric laplace distribution
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
#'If we define the variable Oi, which takes value equal to 1 when ith observation
#'is an outlier, and 0 otherwise, then we propose to calculate the probability of
#'an observation being an outlier as:
#'
#'\deqn{P(O_{i} = 1) = \frac{1}{n-1}\sum{P(v_{i}>v_{j}|data)} \quad (1)}
#'We believe that for points, which are not outliers, this probability should be
#'small, possibly close to zero. Given the natrual ordering of the residuals, it is
#'expected that some observations present greater values for this probability in
#'comparison to others. What we think that should be deemed as an outlier, ought to
#'be those observations with a higher \eqn{P(O_{i} = 1)}, and possibly one that is
#'particularly distant from the others.
#'
#'The probability in the equation can be approximated given the MCMC draws, as follows
#'
#'\deqn{P(O_{i}=1)=\frac{1}{M}\sum{I(v^{(l)}_{i}>max v^{k}_{j})}}
#'
#'where \eqn{M} is the size of the chain of \eqn{v_{i}} after the burn-in period and
#'\eqn{v^{(l)}_{j}} is the \eqn{l}th draw of chain.
#'
#'@references
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'@seealso \code{bayesKL}

bayesProb <- function(y, x, tau, M){
   coefs <- bayesQR::bayesQR(y ~ x, quantile = tau, alasso = FALSE,
                    ndraw = M, prior = NULL)
   beta <- summary(coefs)[[1]]$betadraw[, 1]
   sigma <- summary(coefs)[[1]]$sigmadraw[, 1]
   taup2 <- (2/(tau * (1 - tau)))
   theta <- (1 - 2 * tau) / (tau * (1 - tau))
   x <- cbind(1, x)
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
