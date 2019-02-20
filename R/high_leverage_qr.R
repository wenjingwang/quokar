#' @title High leverage diagnosis for quantile regression model
#' 
#' @description This function provides a better method to detect high leverage observations
#' in quantile regression comparing to the robsut distance method provide in function
#' `frame_distance`. 
#' 
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms,
#' separated by + operators, on the right.
#' @param data data.frame, data used to fit quantile regression model
#' @param tau singular or vector, the quantile(s) to be estimated, strictly between 0 and 1.
#' @return Value of high leverage diagnostic statistics on each quantile; Cut-off value of the 
#' diagnostic statistics; Evaluation of the leverage problem in regression models; Potential 
#' high leverage observation flagged for regression model on each quantile.
#' 
#' @details The proposed statistic use the concept of elemental sets (ESs). ESs consist of 
#' exactly the minimum number of observations to fit a regression mdoel. In linear model setting, 
#' for a sample size $n$ with $p$ predictors, the total number of ERs is $K=C_{n}^{p}$. The $K$ ERs 
#' consist of the set of all feasible solutions to the Linear Programming (LP) problem giving 
#' $RQs$ as optimal solutions. Each RQ is associated with an elemental regression weight. Hence, 
#' A plausible way to deal with RQ hgih leverage diagnosis is via ER. The K ERs consist of the set 
#' of all feasible solutions to the Linear Programming (LP) problem giving RQs as its optimal 
#' solutions. At the optimal solution the $p$ basic observations at level zero correspond to a 
#' specific RQ. 
#' @references 
#' Ranganai E, Van Vuuren J O, De Wet T. Multiple case high leverage diagnosis in regression 
#' quantiles. \ref{Communications in Statistics-Theory and Methods, 2014, 43(16): 3343-3370}.
#' @export

high_leverage_qr <- function(formula, data, tau){
  x <- as.matrix(cbind(1, data[, all.vars(update(formula, 0~.))]))
  y <- data[, all.vars(update(formula, .~ 0))]
  n <- length(y)
  p <- ncol(x)
  ntau <- length(tau)
  k <- 1
  T_J <- rep(0, length(tau))
  h_iJ <- matrix(0, nrow = n - p, ncol = ntau)
  for(i in tau){
    rqs <- rq(formula, data = data, tau = i,  method = "br")
    es <- levels(as.factor((frame_br(rqs, i)$fitting_point$index)))
    x_J <- x[as.numeric(es), ]
    x_I <- x[-as.numeric(es), ]
    H_I <- x_I %*% solve(t(x) %*% x) %*% t(x_I)
    w_J <- det(diag(n - p) - H_I)
    H_IJ <- x_I %*% solve(t(x_J) %*% x_J) %*% t(x_I)
    h_iJ[, k] <- diag(H_IJ)
    T_J[k] <- (w_J * sum(diag(H_IJ)))/(n - p)
    k <- k + 1
  }
  return(leverage_obs = h_iJ)
}


