#' @title Outlier diagnostic for quantile regression model using generalized Cook's distance
#'
#' \code{GCDqr} returns the generalized Cook's distance for each observation deleted 
#' from the quantile regression model.
#' 
#' @description Generalized Cook's distance has been widely used to detect outlying 
#' observations in regression analysis. The idea is to measure the influence
#' of a potential outlier based on the difference between estimation of parameter
#' vector from the whole data set and data set with each observation removed.
#' 
#' @param formula Quantile regression model contains dependent variable and 
#' independent variables 
#' @param data Dataframe. Data sets contains the data of dependent variable and 
#' independent variables
#' @param tau Singular between 0 and 1. Quantiles
#' @param beta Estimation of parameter beta of asymmetric Laplace distribution
#' @param sigma Estimation of parameter sigma of asymmetric Laplace distribuion
#' @return Generalized Cook's distance with each observation deleted from the
#' orginal quantile regression model.
#' 
#'
#' @references
#' Benites L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#' Diagnostics for Quantile Regression Using the Asymmetric Laplace
#' Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#' @export

GCDqr <- function(formula, data, tau, beta, sigma){
  x <- cbind(1, data[, all.vars(update(formula, 0~.))])
  y <- data[, all.vars(update(formula, .~0))]
  n <- length(y)
  p <- ncol(x)
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta)^2/(taup2 * sigma)
  gamma2 <- (2 + thep^2/taup2)/sigma
  muc <- y - x %*% beta
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/
    (besselK(sqrt(delta2* gamma2), 0.5))*(sqrt(delta2 / gamma2))^(-1)
  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/
    (besselK(sqrt(delta2*gamma2), 0.5)) * (sqrt(delta2 / gamma2))
  E1 <- matrix(0, nrow = p, ncol = n)
  for(i in 1:n){
    suma2 <- x[-i,] * c(vchpN[-i] * (y[-i] - x[-i,] %*% beta) - thep)
    E1[,i] <- apply(suma2, 2, sum)/(taup2)
  }
  E2 <- 1: n %>%
    map(function(i) {
      muc_i <- y[-i] - x[-i, ] %*% beta
      sum(3*sigma - (vchpN[-i] * muc_i^2 -
                          2 * muc_i * thep + vchp1[-i] *(thep^2 + 2 * taup2))/taup2)
    })
  E2 <- simplify2array(E2)
  Q1_beta <- E1/sigma
  Q1_sigma <- -E2/(2*sigma^2)
  xM <- c(sqrt(vchpN)) * x
  suma1 <- t(xM) %*% (xM)
  Q2_beta <- -(suma1)/(sigma * taup2)
  Q2_sigma <- 3/(2 * sigma^2) - sum((vchpN * muc^2-2 * muc * thep +
                                        vchp1 * (thep^2 + 2 * taup2)))/(sigma^3*taup2)
  GCD_beta <- 1:n %>%
    map(function(i){
      c(Q1_beta[,i]) %*% solve(-Q2_beta) %*% matrix(Q1_beta[,i], ncol = 1)
    })
  GCD_beta <- simplify2array(GCD_beta)
  GCD_sigma <- 1:n %>%
    map(function(i) {
      Q1_sigma[i] * solve(-Q2_sigma) * Q1_sigma[i]
    })
  GCD_sigma <- simplify2array(GCD_sigma)
  GCD <- GCD_beta + GCD_sigma
  return(GCD)
}










