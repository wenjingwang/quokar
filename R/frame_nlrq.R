#' @title Visualization of fitting process of non-linear quantile regression:
#' interior point algorithm
#' @description This function explore the fitting process of nonlinear
#' quantile regression
#' @param formula non-linear quantile regression model
#' @param data data frame
#' @param tau quantiles
#' @param start the initial value of all parameters to estimate, must be a list
#' @return Weighted observations in non-linear quantile regression model fitting using
#' interior algorithm
#' @author Wenjing Wang \email{wenjingwangr@gmail.com}
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
#' \deqn{f(t*)^{'}d = \sum{|f_i(t*)|}}
#' where \eqn{f(t)=(f_i(t))} and
#' \eqn{J(t)=(\partial f_i(t)/\partial t_j)}.
#' Thus, as proposed by Osborne and Watson(1971), one approach to
#' solving  \eqn{min_{t\in R^{p}} \sum{|f_{i}(t)|}} is to solve a
#' succession of linearized \eqn{l_1} problems minimizing
#' \deqn{\sum |f_{i}(t)-J_{i}(t)^{'}\delta|}
#' @export
#' @examples
#' library(tidyr)
#' library(ggplot2)
#' library(purrr)
#' x <- rep(1:25, 20)
#' y <- SSlogis(x, 10, 12, 2) * rnorm(500, 1, 0.1)
#' Dat <- data.frame(x = x, y = y)
#' formula <- y ~ SSlogis(x, Aysm, mid, scal)
#' nlrq_m <- frame_nlrq(formula, data = Dat, tau = c(0.1, 0.5, 0.9))
#' weights <- nlrq_m$weights
#' m <- data.frame(Dat, weights)
#' m_f <- m %>% gather(tau_flag, value, -x, -y)
#' ggplot(m_f, aes(x = x, y = y)) +
#'   geom_point(aes(size = value, colour = tau_flag)) +
#'   facet_wrap(~tau_flag)

frame_nlrq <- function(formula, data, tau, start){
  ntau <- length(tau)
  n <- nrow(data)
  D_s <- matrix(0, nrow = n, ncol = ntau)
  resid <- matrix(0, nrow = n, ncol = ntau)
  for(i in 1:ntau){
   model <- nlrq_m(formula, data = data, tau = tau[i], trace = FALSE, start)
   D <- model$m$D
   ##turn list into matrix
   D_m <- simplify2array(D)^2
   D_s[,i] <- apply(D_m, 1, mean)
   D_s[,i] <- D_s[, i]/sum(D_s[,i])
   resid[,i] <- model$m$resid()[1:n]
  }
  colnames(D_s) <- paste("tau", tau, sep="")
  colnames(resid) <- paste("tau", tau, sep = "")
  return(list(weights = D_s, resid = resid))
}








