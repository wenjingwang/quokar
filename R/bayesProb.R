#'Mean posterior probability for each observation in
#'Baysian quantile regression model
#'
#'@param y vector, dependent variable in quantile regression
#'
#'@param x matrix, design matrix in quantile regression
#'
#'@param tau quantile
#'
#'@param M MCMC draws
#'
#'@param burn burned MCMC draws
#'
#'@details
#'If we define the variable O_{i}, which takes value equal to 1 when ith observation
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
#'More details please refer to the paper in references
#'
#'@references
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'@seealso \code{bayesKL}
#'@examples
#'\dontrun{
#'ais_female <- subset(ais, Sex == 1)
#'y <- ais_female$BMI
#'x <- ais_female$LBM
#'tau <- 0.5
#'M <- 5000
#'burn <- 1000
#'prob <- bayesProb(y, x, tau, M, burn)
#'case <-  1:100
#'dat <- data.frame(case, prob)
#'ggplot(dat, aes(case, prob))+
#'  geom_point() +
#'  geom_text(data = subset(dat, prob > mean(prob) + 2*sd(prob)),
#‘            aes(label = case), vjust = 0, hjust = 0)
#'}
#'
#'
#'

bayesProb <- function(y, x, tau, M, burn){
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
  res <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      res[i,j] <- sum(v[i, ] > max(v[j, ])) / t
    }
  }
  prob <- rowSums(res) / (n-1)
  return(prob)
}
