#' @title Function to compute minimum covariance determinant and draw
#' mcd-residual plot to diagnose outliers
#' @param object quantile regression model using 'br' or 'fn' method
#' @param tau quantile regression
#'
#' @details The generalized MCD algorithm based on the fast-MCD
#' algorithm formulated by Rousseeuw and Van Driessen(1999), which
#' is similar to the algorithm for least trimmed squares(LTS).
#' The canonical Mahalanobis distance is defined as
#'
#' \deqn{MD(x_i)=[(x_i-\bar{x})^{T}\bar{C}(X)^{-1}(x_i-\bar{x})]^{1/2}}
#' where \eqn{\bar{x}=\frac{1}{n}\sum_{i=1}^{n}x_i} and
#' \eqn{\bar{C}(X)=\frac{1}{n-1}\sum_{i=1}^{n}(x_i-\bar{x})^{T}(x_i-
#' \bar{x})} are the empirical multivariate location and scatter,
#' respectively. Here \eqn{x_i=(x_{i1},...,x_{ip})^{T}} exclueds the
#' intercept. The relation between the Mahalanobis distance
#' \eqn{MD(x_i)} and the hat matrix
#' \eqn{H=(h_{ij})=X(X^{T}X)^{-1}X^{T}} is
#'
#' \deqn{h_{ii}=\frac{1}{n-1}MD^{2}_{i}+\frac{1}{n}}
#'
#' The canonical robust distance is defined as
#'
#' \deqn{RD(x_{i})=[(x_{i}-T(X))^{T}C(X)^{-1}(x_{i}-T(X))]^{1/2}}
#'
#' where \eqn{T(X)} and \eqn{C(X)} are the robust multivariate
#' location and scatter, respectively, obtained by MCD.
#'
#' To achieve robustness, the MCD algorithm estimates the covariance
#' of a multivariate data set mainly through as MCD \eqn{h}-point
#' subset of data set. This subset has the smallest sample-covariance
#' determinant among all the possible \eqn{h}-subsets. Accordingly,
#' the breakdown value for the MCD algorithm equals
#' \eqn{\frac{(n-h)}{n}}. This means the MCD estimates is reliable,
#' even if up to \eqn{\frac{100(n-h)}{n}}% observations in the data
#' set are contaminated.
#' It is possible that the original data is in \eqn{p} dimensional
#' space, but the \eqn{h}-point subset that yields the minimum
#' covariance determinant lies in a lower-dimensional hyperplane.
#' Applying the canonical MCD algorithm to such a data set would result
#' in a singular covariance problem, so that the relevant robust
#' distances cannot be computed. To deal with the singularity problem
#' and provide further leverage point analysis. We implement a
#' generalized MCD algorithm to distinguish in-(hyper) plane points
#' from off-(byper) plane points, and performs MCD leverage point
#' analysis in the dimension-reduced space by projecting all point
#' onto the hyperplane.
#' leverage point and outlier detection
#' the regular variable leverage is defined as
#' \deqn{leverage=\left\{begin{array}{rcl} 0 & \mbox{if RD(x_i) \leq
#' C(p) \\ 1 & \mbox{otherwise}} \end{array}\right}
#' where \eqn{C(p)=\sqrt{\chisq^{2}_{p;1-\alpha}}} is the cutoff value.
#'
#' @description the standardized residuals from quantile regression
#' against the robust MCD distance. This display is used to diagnose
#' both vertical outlier and horizontal leverage points.
#' @author Wenjing Wang
#' @examples
#' library(quantreg)
#' library(ggplot2)
#' library(ALDqr)
#' library(purrr)
#' library(robustbase)
#' library(tidyr)
#' data(ais)
#' tau = c(0.1, 0.5, 0.9)
#' object <- rq(BMI ~ LBM + Ht, data = ais, tau = tau)
#' plot_distance <- frame_distance(object, tau = c(0.1, 0.5, 0.9))
#' distance <- plot_distance[[1]]
#' cutoff_v <- plot_distance[[2]]
#' cutoff_h <- plot_distance[[3]]
#' n <- nrow(object$model)
#' case <- rep(1:n, length(tau))
#' distance <- cbind(case, distance)
#' ggplot(distance, aes(x = rd, y = value)) +
#'    geom_point() +
#'    geom_vline(xintercept = cutoff_v, colour = "red") +
#'    geom_hline(yintercept = cutoff_h, colour = "red") +
#'    facet_wrap(~ tau_flag, scale = 'free_x') +
#'    geom_text(data = subset(distance, value > cutoff_h[1] |
#'                                      value < cutoff_h[2] |
#'                                      rd > cutoff_v),
#'              aes(label = case)) +
#'    xlab("Robust Distance") +
#'    ylab("Residuals")

#' ggplot(distance, aes(x = md, y = value)) +
#'    geom_point() +
#'    geom_vline(xintercept = cutoff_v, colour = "red") +
#'    geom_hline(yintercept = cutoff_h, colour = "red") +
#'    facet_wrap(~ tau_flag, scale = 'free_x') +
#'    geom_text(data = subset(distance, value > cutoff_h[1] |
#'                                      value < cutoff_h[2] |
#'                                      rd > cutoff_v),
#'              aes(label = case)) +
#'    xlab("Mahalanobis Distance") +
#'    ylab("Residuals")

frame_distance <- function(object, tau){
    x <- object$model[, -1]
    p <- ncol(x)
    n <- nrow(x)
    resid <- object$residuals
    center_m <- apply(x, 2, mean)
    cov_m <- cov(x)
    md <- rep(0, n)
    for(i in 1:n){
        mid_m <- x[i,] - center_m
        mid_m <- as.matrix(mid_m)
        md[i] <- sqrt(mid_m%*%solve(cov_m)%*%matrix(mid_m,
                                                    ncol= 1))
    }
    md <- matrix(md, ncol = 1)
    mcd <- covMcd(x)
    Tx <- mcd$center
    Cx <- mcd$cov
    rd <- rep(0, n)
    for(i in 1:n){
        mid_r <- x[i, ] - Tx
        mid_r <- as.matrix(mid_r)
       rd[i] <- sqrt(mid_r%*%solve(Cx)%*%matrix(mid_r,ncol = 1))
    }
    rd <- matrix(rd, ncol = 1)
    cutoff_v <- qchisq(p = 0.975, df = p)
    cutoff_h <- c(-2.5, 2.5)
    ##data set for plotting
    rd_m <- data.frame(md, rd, resid)
    rd_f <- rd_m %>% gather(tau_flag, value, -md, -rd)
    return(list("Distance" = rd_f, "Vertical Cutoff" = cutoff_v,
           "Horizental Cutoff" = cutoff_h))
}















