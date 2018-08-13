#' @title This function provide visualization method for the quantile 
#' regression model.
#' 
#' \code{qrdepth.plot} returns the visualization of outliers in quantile
#' regression model estimating.
#' 
#' @description This function provides visualization of outliers in quantile regression
#' depend on depth estimator.
#' 
#' @references Rousseeuw P J, Hubert M. Regression depth. 
#' \emph{Journal of the American Statistical Association, 1999, 94(446): 388-402}.
#' 
#' Debruyne M, Hubert M, Portnoy S, et al. Censored depth quantiles. 
#' \emph{Computational statistics & data analysis, 2008, 52(3): 1604-1614}.
#' 
#' @examples
#' library(ggplot2)
#' data("CYGOB1", package = "HSAUR")
#' head(CYGOB1)
#' ggplot(CYGOB1, aes(logst, logli)) +
#'   geom_point() +
#'   geom_quantile(quantiles = seq(0.1:0.9, by = 0.1))
#' ## depth estimation
#' formula <- logli ~ logst
#' data <- CYGOB1
#' qrdepth1 <- qrdepth(formula, data, tau = 0.1)
#' qrdepth2 <- qrdepth(formula, data, tau = 0.2)
#' qrdepth3 <- qrdepth(formula, data, tau = 0.3)
#' qrdepth4 <- qrdepth(formula, data, tau = 0.4)
#' qrdepth5 <- qrdepth(formula, data, tau = 0.5)
#' qrdepth6 <- qrdepth(formula, data, tau = 0.6)
#' qrdepth7 <- qrdepth(formula, data, tau = 0.7)
#' qrdepth8 <- qrdepth(formula, data, tau = 0.8)
#' qrdepth9 <- qrdepth(formula, data, tau = 0.9)
#' beta1 <- qrdepth1$beta
#' beta2 <- qrdepth2$beta
#' beta3 <- qrdepth3$beta
#' beta4 <- qrdepth4$beta
#' beta5 <- qrdepth5$beta
#' beta6 <- qrdepth6$beta
#' beta7 <- qrdepth7$beta
#' beta8 <- qrdepth8$beta
#' beta9 <- qrdepth9$beta
#' betas <- data.frame(a = c(beta1[2], beta2[2], beta3[2],
#'                           beta4[2], beta5[2], beta6[2],
#'                           beta7[2], beta8[2], beta9[2]),
#'                     b = c(beta1[1], beta2[1], beta3[1],
#'                           beta4[1], beta5[1], beta6[1],
#'                           beta7[1], beta8[1], beta9[1]))
#' plot.outlier.qrdepth(qrdepth1)


plot.outlier.qrdepth <- function(qrdepth.object){
  if(class(qrdepth.object) != "qrdepth"){
    stop(simpleError("The argument 'qrdepth.object' is not an object 
                     of class 'qrdepth'"))
  } 
  x <- qrdepth.object$x
  y <- qrdepth.object$y
  tau <- qrdepth.object$tau
  beta <- qrdepth.object$beta
  residuals <- y - as.matrix(x) %*% as.matrix(qrdepth.object$beta)
  case <- 1:length(y)
  resid_dat <- data.frame(case, residuals)
  ggplot(resid_dat, aes(x = case, y = residuals)) +
    geom_point() +
    geom_text(data = dplyr::filter(resid_dat, residuals > mean(residuals) +
                                     2*sd(residuals)),
              aes(label = case), vjust = 1)+
    theme(aspect.ratio = 1)
}












