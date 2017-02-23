#' @title Weighting Matrix of Quantile regression using Asymmetric Laplace Distrubtion
#' @description This function calulate the weighting matrix
#' @param y dependent variable of quantile regression
#' @param x design matrix of quantile regression
#' @param tau quantile must be a scaler
#' @param error The EM algorithm accuracy of error used in MLE estimation
#' @param iter The iteration frequancy for EM algorithm used in MLE estimation
#'
#' @details In the estimation procedure in EM algorithm, we can see that
#' \eqn{\varepsilon} is inversely proportional to
#' \eqn{d_i = |y_i-x^{'}_{i}\beta^{(k)}_{p}|/\sigma}.
#' Hence, \eqn{u_i(\theta^{k})=\varepsilon_{-1i}(\theta^{(k)})}
#' can be interpreted as a type of weight for \eqn{i}th case in the estimates of
#' \eqn{\beta_{(k)^p}}, which tends to be small for outlying observations.
#'
#' @author Wenjing Wang
#' @export
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' library(ALDqr)
#' data(ais)
#' y <- ais$BMI
#' x <- cbind(1, ais$LBM)
#' tau <-  c(0.1, 0.5, 0.9)
#' error <- 1e-06
#' iter <- 100
#' weights <- frame_ald_weight(y, x, tau, error, iter)
#' weights
#'
frame_ald_weight <- function(y, x, tau, error, iter){
    ntau <- length(tau)
    n <- length(y)
    p <- ncol(x)
    vchpN <- matrix(0, nrow = n, ncol = ntau)
    for(i in 1:ntau){
    qr <- EM.qr(y,x,tau[i],error,iter)
    beta_qr <- qr$theta[1:p,]
    sigma_qr <- qr$theta[p+1]
    taup2 <- (2/(tau[i] * (1 - tau[i])))
    thep <- (1 - 2 * tau[i]) / (tau[i] * (1 - tau[i]))
    delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
    gamma2 <- (2 + thep^2/taup2)/sigma_qr
    vchpN[, i] <- (besselK(sqrt(delta2 * gamma2), 0.5 - 1)/
      (besselK(sqrt(delta2 * gamma2), 0.5))) *
      (sqrt(delta2 / gamma2))^(-1)
    vchpN[, i] <- vchpN[,i]/sum(vchpN[,i])
    }
    colnames(vchpN) <- paste("tau", tau, sep = "")
    return(vchpN)
}

