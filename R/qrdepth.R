#' @title Regression depth estimator for quantile regression model
#'
#' \code{qrdepth} returns the quantile regression depth estimator for quantile
#' regression model
#' 
#' @description Regression depth estimator is more robsut comparing to the mle 
#' estimator. The geometrical interpretation of this estimator is the 
## smallest number of observations one has to pass in order to turn
## the hyperplane $\beta$ into vertical position.
#' 
#' @param formula Quantile regression model contains dependent variable and 
#' independent variables 
#' @param data Dataframe. Data sets contains the data of dependent variable and 
#' independent variables
#' @param tau Singular between 0 and 1. Quantiles
#' @param Initial Initial size of combinations of observations for estimating 
#' quantile regression model from $L_{1}$ algorithm
#' @param Iteration Replacing size of combinations observations for estimating 
#' quantile regression model from $L_{1}$ algorithm
#' @return Regression depth estimator for quantile regression model
#' @references Rousseeuw P J, Hubert M. Regression depth. 
#' \emph{Journal of the American Statistical Association, 1999, 94(446): 388-402}.
#' 
#' Debruyne M, Hubert M, Portnoy S, et al. Censored depth quantiles. 
#' \emph{Computational statistics & data analysis, 2008, 52(3): 1604-1614}.
#' 
#' @export
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
#' beta1 <- qrdepth(formula, data, tau = 0.1)$beta
#' beta2 <- qrdepth(formula, data, tau = 0.2)$beta
#' beta3 <- qrdepth(formula, data, tau = 0.3)$beta
#' beta4 <- qrdepth(formula, data, tau = 0.4)$beta
#' beta5 <- qrdepth(formula, data, tau = 0.5)$beta
#' beta6 <- qrdepth(formula, data, tau = 0.6)$beta
#' beta7 <- qrdepth(formula, data, tau = 0.7)$beta
#' beta8 <- qrdepth(formula, data, tau = 0.8)$beta
#' beta9 <- qrdepth(formula, data, tau = 0.9)$beta
#' betas <- data.frame(a = c(beta1[2], beta2[2], beta3[2],
#'                           beta4[2], beta5[2], beta6[2],
#'                           beta7[2], beta8[2], beta9[2]),
#'                     b = c(beta1[1], beta2[1], beta3[1],
#'                           beta4[1], beta5[1], beta6[1],
#'                           beta7[1], beta8[1], beta9[1]))
#' ggplot(CYGOB1, aes(logst, logli)) +
#'    geom_point() +
#'    geom_quantile(quantiles = seq(0.1, 0.9, by = 0.1)) +
#'    geom_abline(slope = betas$a, intercept = betas$b, colour = "red")

qrdepth <- function(formula, data, tau, Initial = 10000, Iteration = 10000){
  x <- cbind(1, data[, all.vars(update(formula, 0~.))])
  y <- data[, all.vars(update(formula, .~ 0))]
  n <- length(y)
  p <- ncol(x)
  M <- Initial + Iteration
  betaseries <- matrix(0, nrow = p, ncol = M)
  for(i in 1:M){
    k <- 0
    while(k == 0){
      a <- sample(1:n, p)
      x_square <- x[a, ]
      if(qr(x_square)$rank == p) k <- 1
    }
    betai <- solve(x_square) %*% y[a]
    betaseries[, i] <- betai
  }
  Initial_beta <- betaseries[, 1:Initial]
  Iteration_beta <- betaseries[, (Initial + 1): M]
  beta <- maxObject(x, y, tau, Initial_beta, Iteration_beta)$beta
  ret <-  list()
  ret <- list(beta = beta, tau = tau, x = x, y = y)
  class(ret) <-  "qrdepth"
  return(ret)
}

maxObject <- function(x, y, tau, Initial_beta, Iteration_beta){
  x <- as.matrix(x)
  n <- length(y)
  Initial_num <- ncol(Initial_beta)
  Iteration_num <- ncol(Iteration_beta)
  nc <- 1:n
  i <- 1
  maxobjf <- 0
  while(i <= Initial_num){
    betai <- Initial_beta[, i]
    ri <- y - as.matrix(x) %*% matrix(betai, nrow = p, ncol = 1)
    minobjBetai <- 1e+20
    j <- 1
    while(j <= Initial_num & minobjBetai > maxobjf){
      Iteration_betaj <- betai - Iteration_beta[, j]
      if(length(which(abs(Iteration_betaj) > 1e-07)) != 0){
        proj <- as.matrix(x) %*% Iteration_betaj
        Kcomp <- 1:n
        riKcompos <- ri[Kcomp] >= -1e-07
        projKcompos <- proj[Kcomp] > 1e-07
        riKcompne <- ri[Kcomp] <= 1e-07
        projKcompne <- proj[Kcomp] < -1e-07
        obj1V <- tau * sum(riKcompos * projKcompne) +
          (1 - tau) * sum(riKcompne * projKcompos)
        obj1G <- tau * sum(riKcompos * projKcompos) +
          (1 - tau) * sum(riKcompne * projKcompne)
        minobjBetai <- min(obj1G, obj1V, minobjBetai)
      }
      j <- j + 1
    }
    if( minobjBetai > maxobjf){
      beta <- betai
      maxobjf <- minobjBetai
    }
    i <- i + 1
  }
  return(list(beta = beta, maxobjf = maxobjf))
}

