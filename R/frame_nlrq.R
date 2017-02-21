#' @title Explore fitting process of non-linear quantile  regression
#' @description
#' @param object non-linear quantile regression model
#' @param data data frame
#' @param tau quantiles
#' @author Wenjing Wang
#' @details To extentd the linear programming method to the case of
#' non-linear response functions, Koenker & Park(1996) considered the
#' nonlinear \eqn{l_{1}} problem
#' \deqn{min_{t\in R^{p}} \sum{|f_{i}(t)|}}
#' where, for example,
#' \deqn{f_{i}(t)=y_i-f_{0}(x_i, t)}
#' As noted by El Attar et al(1979) a necessary condition for \eqn{t*}
#' to solve \eqn{min_{t\in R^{p}} \sum{|f_{i}(t)|}} is that there
#' exists a vector \eqn{d \in [-1, 1]^n} such that
#' \deqn{J(t*)^{'}d = 0}
#' \deqn{f(t*)^{'}d = \sum{|f_i(t*)|}
#' where \eqn{f(t)=(f_i(t))} and
#' \eqn{J(t)=(\partial f_i(t)/\partial t_j)}.
#' Thus, as proposed by Osborne and Watson(1971), one approach to
#' solving  \eqn{min_{t\in R^{p}} \sum{|f_{i}(t)|}} is to solve a
#' succession of linearized \eqn{l_1} problems minimizing
#' \deqn{\sum |f_{i}(t)-J_{i}(t)^{'}\delta|}
#' @export
#'
frame_nlrq <- function(object, data, tau){
   D2 <- quokar::nlrq(formula, data = data , tau = tau,trace = FALSE,
         method = L-BFGS-B)
   resid <- object$residuals
}




