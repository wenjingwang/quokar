globalVariables("tau_flag")
#' @title Residual-robust distance plot of quantile regression model
#' @description the standardized residuals from quantile regression
#' against the robust MCD distance. This display is used to diagnose
#' both vertical outlier and horizontal leverage points. Function
#' `frame_distance` only work for linear quantile regression model. With
#' non-linear model, use `frame_distance_complex`
#' @param x matrix, covariate of quantile regression model
#' @param resid matrix, residuals of quantile regression models
#' @param tau singular or vectors, quantile
#' @return dataframe for residual-robust distance plot
#' @details The generalized MCD algorithm based on the fast-MCD
#' algorithm formulated by Rousseeuw and Van Driessen(1999), which
#' is similar to the algorithm for least trimmed squares(LTS).
#' The canonical Mahalanobis distance is defined as
#' \deqn{MD(x_i)=[(x_i-\bar{x})^{T}\bar{C}(X)^{-1}(x_i-\bar{x})]^{1/2}}
#' where \eqn{\bar{x}=\frac{1}{n}\sum_{i=1}^{n}x_i} and
#' \eqn{\bar{C}(X)=\frac{1}{n-1}\sum_{i=1}^{n}(x_i-\bar{x})^{T}(x_i-
#' \bar{x})} are the empirical multivariate location and scatter,
#' respectively. Here \eqn{x_i=(x_{i1},...,x_{ip})^{T}} exclueds the
#' intercept. The relation between the Mahalanobis distance
#' \eqn{MD(x_i)} and the hat matrix
#' \eqn{H=(h_{ij})=X(X^{T}X)^{-1}X^{T}} is
#' \deqn{h_{ii}=\frac{1}{n-1}MD^{2}_{i}+\frac{1}{n}}
#' The canonical robust distance is defined as
#' \deqn{RD(x_{i})=[(x_{i}-T(X))^{T}C(X)^{-1}(x_{i}-T(X))]^{1/2}}
#' where \eqn{T(X)} and \eqn{C(X)} are the robust multivariate
#' location and scatter, respectively, obtained by MCD.
#' To achieve robustness, the MCD algorithm estimates the covariance
#' of a multivariate data set mainly through as MCD \eqn{h}-point
#' subset of data set. This subset has the smallest sample-covariance
#' determinant among all the possible \eqn{h}-subsets. Accordingly,
#' the breakdown value for the MCD algorithm equals
#' \eqn{\frac{(n-h)}{n}}. This means the MCD estimates is reliable,
#' even if up to \eqn{\frac{100(n-h)}{n}}% observations in the data
#' set are contaminated.
#' @author Wenjing Wang \email{wenjingwangr@gmail.com}
#' @export
frame_distance_complex <- function(x,resid,tau){
  p <- ncol(x)
  n <- nrow(x)
  colnames(resid) <- paste("tau", tau, sep="")
  center_m <- apply(x, 2, mean)
  cov_m <- cov(x)
  md <- rep(0, n)
  for(i in 1:n){
    mid_m <- x[i,] - center_m
    mid_m <- as.matrix(mid_m)
    md[i] <- sqrt(matrix(mid_m, nrow = 1)%*%solve(cov_m)%*%matrix(mid_m, ncol= 1))
  }
  md <- matrix(md, ncol = 1)
  mcd <- covMcd(x)
  Tx <- mcd$center
  Cx <- mcd$cov
  rd <- rep(0, n)
  for(i in 1:n){
    mid_r <- x[i, ] - Tx
    mid_r <- as.matrix(mid_r)
    rd[i] <- sqrt(matrix(mid_r, nrow = 1)%*%solve(Cx)%*%matrix(mid_r,ncol = 1))
  }
  rd <- matrix(rd, ncol = 1)
  cutoff_v <- sqrt(qchisq(p = 0.975, df = p))
  cutoff_h <- rep(0, length(tau))
  for(i in 1:length(tau)){
    cutoff_h[i] <- median(abs(resid[,i])/qnorm(p=0.75, mean=0, sd = 1))
  }
  cutoff_h <- 5*cutoff_h
  rd_m <- data.frame(md, rd, resid)
  rd_f <- rd_m %>% gather(tau_flag, residuals, -md, -rd)
  return(list("Distance" = rd_f, "Vertical Cutoff" = cutoff_v,
              "Horizental Cutoff" = cutoff_h))
}

