#' @title Outlier diagnostic for quantile regression model using Kullback¨CLeibler divergence
#'
#' \code{bayesKL_qr} returns the Kullback¨CLeibler divergence for each observation 
#' in quantile regression model.
#' 
#' @description Kullback-Leibler divergence is a method of measuring the
#' distance between those latent variables distribution in the Baysian quantile 
#' regression framework. 
#' 
#' @param formula Quantile regression model contains dependent variable and 
#' independent variables 
#' @param data Dataframe. Data sets contains the data of dependent variable and 
#' independent variables
#' @param tau Singular between 0 and 1. Quantiles
#' @param beta Estimation of parameter beta of asymmetric Laplace distribution
#' @param sigma Estimation of parameter sigma of asymmetric Laplace distribuion
#' @param ndraw Integer. MCMC draws
#' @return Kullback¨CLeibler divergence of each observation in regression model.
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
#' Outlier diagnostic and visualization
#' ## tau = 0.1
#' kl <- bayesKL_qr(formula = BMI ~ LBM, data = ais_female, tau = 0.1,
#'              beta = beta, sigma = sigma, ndraw = 500)
#' kl_dat <- data.frame(case = 1:nrow(ais_female), kl = kl)
#' p1 <- ggplot(kl_dat, aes(x = case, y = kl)) +
#'   geom_point() +
#'   geom_point(data = subset(kl_dat, kl > median(kl) + 3 * sd(kl)), 
#'             colour = "red") +
#'   geom_text(data = subset(kl_dat, kl > median(kl) + 3 * sd(kl)), 
#'             aes(label = case))
#' ## tau = 0.5
#' kl <- bayesKL_qr(formula = BMI ~ LBM, data = ais_female, tau = 0.5,
#'               beta = beta, sigma = sigma, ndraw = 500)
#' kl_dat <- data.frame(case = 1:nrow(ais_female), kl = kl)
#' p2 <- ggplot(kl_dat, aes(x = case, y = kl)) +
#'   geom_point() +
#'   geom_point(data = subset(kl_dat, kl > median(kl) + 3 * sd(kl)), 
#'              colour = "red") +
#'   geom_text(data = subset(kl_dat, kl > median(kl) + 3 * sd(kl)), 
#'              aes(label = case))
#' ## tau = 0.9
#' kl <- bayesKL_qr(formula = BMI ~ LBM, data = ais_female, tau = 0.9,
#'               beta = beta, sigma = sigma, ndraw = 500)
#' kl_dat <- data.frame(case = 1:nrow(ais_female), kl = kl)
#' p3 <- ggplot(kl_dat, aes(x = case, y = kl)) +
#'   geom_point() +
#'   geom_point(data = subset(kl_dat, kl > median(kl) + 3 * sd(kl)), 
#'              colour = "red") +
#'   geom_text(data = subset(kl_dat, kl > median(kl) + 3 * sd(kl)), 
#'             aes(label = case))
#' ## visualization
#' grid.arrange(p1, p2, p3, ncol = 3)
#' export



bayesKL_qr <- function(formula, data, tau, beta, sigma, ndraw = 1000){
  x <- as.matrix(cbind(1, data[, all.vars(update(formula, 0~.))]))
  y <- data[, all.vars(update(formula, .~ 0))]
  n <- length(y)
  p <- ncol(x)
  taup2 <- (2/(tau * (1 - tau)))
  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  v <- matrix(0, nrow = n, ncol = ndraw)
  for (i in 1:n) {
    param1 <- 1/2
    param2 <- (y[i] - x[i, ] %*% t(beta))^2/(taup2 * sigma)
    param3 <- 2/sigma + theta^2/(taup2 * sigma)
    v[i, ] <- GVGrvg::rgig(ndraw, param1, param2, param3)
  }
  KLS <- matrix(0, nrow = n, ncol = n)
  hs <- rep(0, n)
  for(i in 1:n){
    hs[i] <- stats::density(v[i,], kernel = "gaussian")$bw
  }
  g <- matrix(0, nrow = n, ncol = ndraw)
  for(i in 1:n){
    hi <- hs[i]
    upper_x <- max(v)
    lower_x <- min(v)
    ranges <- seq(lower_x, upper_x, length.out = ndraw)
    for(j in 1:ndraw){
      g[i, j] <- gaussian_kernel(ranges[j], v[i, ], ndraw, hi)
    }
  }
  for(i in 1:n) {
    for(j in 1:n) {
      internal_cal <- log(g[i, ]/g[j, ]) * g[i, ]
      KLS[i,j] <- trapz(ranges, internal_cal)
    }
  }
  diag(KLS) <- 0
  KLD <- rowSums(KLS)/(n-1)
  return(KLD)
}
gaussian_kernel <- function(x, y, ndraw, h){
  mid <- ((x - y)/h)^2
  fun1 <- sum((1/sqrt(2 * pi)) * exp(-mid/2))/(ndraw * h)
  return(fun1)
}
trapz <- function(x, y){
  idx <- 2:length(x)
  return (as.double((x[idx]-x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
}
