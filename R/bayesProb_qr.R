#' @title Outlier diagnostic for quantile regression model using mean probability posterior
#'
#' \code{bayesProb_qr} returns the mean probability posterior for each observation 
#' in quantile regression model.
#' 
#' @description Mean probability posterior is a method of measuring the
#' probability of being outliers of each observation in quatnile regression
#' model.  
#' 
#' @param formula Quantile regression model contains dependent variable and 
#' independent variables 
#' @param data Dataframe. Data sets contains the data of dependent variable and 
#' independent variables
#' @param tau Singular between 0 and 1. Quantiles
#' @param beta Estimation of parameter beta of asymmetric Laplace distribution
#' @param sigma Estimation of parameter sigma of asymmetric Laplace distribuion
#' @param ndraw Integer. MCMC draws
#' @return Mean probability of each observation in regression model.
#' 
#' @references  Santos B, Bolfarine H. (2016). On Bayesian quantile regression and 
#' outliers.\emph{arXiv preprint arXiv:1601.07344}
#' 
#' @examples
#' library(tidyr)
#' library(dplyr)
#' library(ggplot2)
#' library(bayesQR)
#' library(GIGrvg)
#' library(gridExtra)
#' ## Baysian quantile regression estimator and parameters of asymmetric Laplace
#' ## distribution
#' data(ais)
#' ais_female <- filter(ais, Sex == 1)
#' bayes_qr <- bayesQR(BMI ~ LBM, data = ais_female, quantile = 0.1, 
#'                    alasso = TRUE, ndraw = 500)
#' beta <- matrix(summary(bayes_qr)[[1]]$betadraw[,1])
#' sigma <- summary(bayes_qr)[[1]]$sigmadraw[1]
#' ## Outlier diagnostics and visualization
#' ## tau = 0.1
#' probs <- bayesProb_qr(formula = BMI ~ LBM, data = ais_female, tau = 0.1, 
#'                    beta = beta, sigma = sigma, ndraw = 1000)
#' probs_dat <- data.frame(case = 1:nrow(ais_female), probs = probs)
#' p1 <- ggplot(probs_dat, aes(x = case, y = probs)) +
#'    geom_point() +
#'    geom_point(data = subset(probs_dat, probs > median(probs) + 3 * sd(probs)), 
#'               colour = "red") +
#'    geom_text(data = subset(probs_dat, probs > median(probs) + 3 * sd(probs)), 
#'              aes(label = case))
#' ## tau = 0.5
#' probs <- bayesProb_qr(formula = BMI ~ LBM, data = ais_female, tau = 0.5, 
#'                    beta = beta, sigma = sigma, ndraw = 1000)
#' probs_dat <- data.frame(case = 1:nrow(ais_female), probs = probs)
#' p2 <- ggplot(probs_dat, aes(x = case, y = probs)) +
#'     geom_point() +
#'     geom_point(data = subset(probs_dat, probs > median(probs) + 3 * sd(probs)), 
#'                colour = "red") +
#'     geom_text(data = subset(probs_dat, probs > median(probs) + 3 * sd(probs)), 
#'               aes(label = case))
#' ## tau = 0.9
#' probs <- bayesProb_qr(formula = BMI ~ LBM, data = ais_female, tau = 0.9, 
#'                   beta = beta, sigma = sigma, ndraw = 1000)
#' probs_dat <- data.frame(case = 1:nrow(ais_female), probs = probs)
#' p3 <- ggplot(probs_dat, aes(x = case, y = probs)) +
#'   geom_point() +
#'   geom_point(data = subset(probs_dat, probs > median(probs) + 3 * sd(probs)), 
#'              colour = "red") +
#'   geom_text(data = subset(probs_dat, probs > median(probs) + 3 * sd(probs)), 
#'             aes(label = case))
#' ## visualization
#' grid.arrange(p1, p2, p3, ncol = 3)
#' export




bayesProb_qr <- function(formula, data, tau, beta, sigma, ndraw){
  x <- as.matrix(cbind(1, data[, all.vars(update(formula, 0~.))]))
  y <- data[, all.vars(update(formula, .~ 0))]
  n <- length(y)
  taup2 <- (2/(tau * (1 - tau)))
  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  v <- matrix(0, nrow = n, ncol = ndraw)
  for (i in 1:n) {
      param1 <- 1/2
      param2 <- (y[i] - x[i, ] %*% t(beta))^2/(taup2 * sigma)
      param3 <- 2/sigma + theta^2/(taup2 * sigma)
      v[i, ] <- rgig(ndraw, param1, param2, param3)
  }
  res <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      res[i,j] <- sum(v[i, ] > max(v[j, ])) / ndraw
    }
  }
  prob <- rowSums(res) / (n-1)
  return(prob)
}
