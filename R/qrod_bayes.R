#'@title Outlier Dignostic for Quantile Regression Based on Bayesian Estimation
#'
#'@param y dependent variable in quantile regression
#'
#'@param x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix
#'
#'@param tau quantile
#'
#'@param M the iteration frequancy for MCMC used in Baysian Estimation
#'
#'@param method the diagnostic method for outlier detection
#'
#'@description
#'This function cacluate the mean probability of posterior
#'of quantile regression model with asymmetric laplace distribution based on bayes
#'estimation procedure.
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
#'Another proposal to address these differences between the posterior distributions
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
#'@author Wenjing Wang \email{wenjingwang1990@gmail.com}
#'@keywords quantile regression outlier diagnostics
#'@references Benites L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#' Diagnostics for Quantile Regression Using the Asymmetric Laplace
#' Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#'Hawkins D M, Bradu D, Kass G V.(1984)``Location of several outliers in
#'multiple-regression data using elemental sets. \emph{Technometrics},
#'26(3), 197-208.
#'
#'Koenker R, Bassett Jr G.(1978)`` Regression quantiles,
#'\emph{Econometrica}, 1, 33-50.
#'
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'Kozumi H, Kobayashi G.(2011)``Gibbs sampling methods for Bayesian
#'quantile regression,\emph{Journal of statistical computation and
#'simulation}, 81(11), 1565-1578.
#'
#'@seealso \code{qrod_mle}
#'@examples
#'\dontrun{
#'library(ALDqr)
#'data(ais)
#'y <- ais$BMI
#'sexInd <- (ais$Sex == 1) + 0
#'x <- cbind(ais$LBM, ais$sexInd)
#'qrod_bayes(y, x, tau = 0.1, M = 100, method = "bayes.prob")
#'qrod_bayes(y, x, tau = 0.1, M = 100, method = "bayes.kl")
#'}
#'
#'@export
#'
#'
qrod_bayes <- function(y, x, tau, M,
                       method = c("bayes.prob", "bayes.kl")){
    method <- match.arg(method)
    distance <- list()
    ntau <- length(tau)
  if(!(method %in% c("bayes.prob","bayes.kl"))){
    stop("Method should be 'bayes.prob' or 'bayes.kl'")
  }
  if(method == "bayes.prob"){
      1 : ntau %>%
          map(function(i){
            names(distance[i]) <- paste('tau=', tau[i], sep = '')
            distance[[i]] <- bayesProb(y, x, tau[i], M)
              })
  }else if(method == "bayes.kl"){
      1 : ntau %>%
          map(function(i){
          names(distance[i]) <- paste('tau=', tau[i], sep = '')
          distance[[i]] <- bayesProb(y, x, tau[i], M)
              })
  }
}




