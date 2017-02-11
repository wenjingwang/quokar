#globalVariables(c("obs", "index", "variable", "value", "tcrit", "alpha"))
#' @title Observations used in quantile regression fitting using br algorithem
#'
#' @details This is a function that can be used to create point plot
#' for the observations used in quantile regression fitting based
#' on 'br'method.
#'
#' @param object quantile regression model using br method
#' @param tau quantiles can be a single quantile or a vector of
#' quantiles
#' @return All observations and the observations used in quantile
#' regression fitting using br algorithem
#' @description get the observation used in br algorithem
#' @examples
#' library(ggplot2)
#' library(quantreg)
#' data(ais)
#' tau <- c(0.1, 0.5, 0.9)
#' object1 <- rq(BMI ~ LBM, tau, method = 'br', data = ais)
#' data_plot <- frame_br(object1, tau)$data_plot
#' choose <- frame_br(object1, tau)$choose
#' ggplot(data_plot,
#'  aes(x=value, y=data_plot[,2])) +
#'  geom_point(alpha = 0.1) +
#'  ylab('y') +
#'  xlab('x') +
#'  facet_wrap(~variable, scales = "free_x", ncol = 2) +
#'  geom_point(data = choose, aes(x = value, y = y,
#'                                       group = tau_flag,
#'                                       colour = tau_flag,
#'                                       shape = obs))
#'
#' object2 <- rq(BMI ~ Ht + LBM + Wt, tau, method = 'br',
#'             data = ais)
#' data_plot <- frame_br(object1, tau)$data_plot
#' choose <- frame_br(object1, tau)$choose
#' ggplot(data_plot,
#'  aes(x=value, y=data_plot[,2])) +
#'  geom_point(alpha = 0.1) +
#'  ylab('y') +
#'  xlab('x') +
#'  facet_wrap(~variable, scales = "free_x", ncol = 2) +
#'  geom_point(data = choose, aes(x = value, y = y,
#'                                       group = tau_flag,
#'                                       colour = tau_flag,
#'                                       shape = obs))
#' @export
#' @useDynLib
#'
frame_br <- function(object, tau){
  y <- matrix(object$y, ncol = 1)
  colnames(y) <-'y'
  x <- object$x
  x <- as.matrix(x)
  ntau <- length(tau)
  h <- matrix(0, nrow = ntau, ncol = ncol(x))
  for (i in 1:ntau){
    h[i, ] <- wh(object, tau[i])
  }
  colnames(h) <- paste('indice', 1:ncol(h), sep='')
  tau_flag <- paste('tau', tau, sep = '=')
  h <- cbind(tau_flag, data.frame(h))
  print('Observations used in br method fitting')
  print(h)
  if(colnames(object$x)[1] == '(Intercept)'){
    x <- object$x[,-1]
    x <- as.matrix(x)
  }
  colnames(x) <- paste('x', 1:ncol(x), sep='')
  data_plot <- data.frame(index = 1:length(y), y, x)
  data_plot_g <- data_plot %>% gather(variable, value, -c(1,2))
  choose_point <- h %>% gather(obs, index, -tau_flag)
  merge_x_y <- choose_point %>%
    inner_join(data_plot, by ='index')
  choose_point2 <- merge_x_y %>% gather(variable, value,
                                        -c(index, tau_flag, obs, y))

  return(list(data_plot = data_plot_g, choose = choose_point2))
}
