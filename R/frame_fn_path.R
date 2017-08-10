#globalVariables(c("variable", "value"))
#' @title Visualization of the fitting path
#' of quantile regression: interior point method
#' @param object quantile regression model using interior point method
#' @param tau quantile
#' @return The fitting path of quantile regression model using interior
#' point method
#' @details This function used to illustrate the fitting process of
#' quantile regression using interior point method.
#' Koenker and Bassett(1978) introduced asymmetric weight on positive
#' and negative residuals, and solves the slightly modified l1-problem.
#' @description observations used in quantile regression fitting
#'
#' \deqn{min_{b \in R^{p}}\sum_{i=1}^{n}\rho_{\tau}(y_i-x_{i}^{'}b)}
#'
#' where \eqn{\rho_{\tau}(r)=r[\tau-I(r<0)]$ for $\tau \in (0,1)}. This
#' yields the modified linear program
#'
#' \deqn{min(\tau e^{'}u+(1-\tau)e^{'}v|y=Xb+u-v,(u,v) \in
#' R_{+}^{2n})}
#'
#' Adding slack variables, \eqn{s}, satisfying the constrains
#' \eqn{a+s=e}, we
#' obtain the barrier function
#'
#' \deqn{B(a, s, u) = y^{'}a+\mu \sum_{i=1}^{n}(loga_{i}+logs_{i})}
#'
#' which should be maximized subject to the constrains
#' \eqn{X^{'}a=(1-\tau)X^{'}e} and \eqn{a+s=e}. The Newton step
#' \eqn{\delta_{a}}
#' solving
#'
#' \deqn{max{y^{'}\delta_{a}+\mu \delta^{'}_{a}(A^{-1}-S^{-1})e-
#' \frac{1}{2}\mu \delta_{a}^{'}(A^{-2}+S^{-2})\delta_{a}}}
#'
#' subject to \eqn{X{'}\delta_{a}=0}, satisfies
#'
#' \deqn{y+\mu(A^{-1}-S^{-1})e-\mu(A^{-2}+S^{-2})\delta_{a}=Xb}
#'
#' for some \eqn{b\in R^{p}}, and \eqn{\delta_{a}} such that
#' \eqn{X^{'}\delta_{a}=0}.
#' Using the constraint, we can solve explicitly for the vector
#' \eqn{b},
#'
#' \deqn{b=(X^{'}WX)^{-1}X^{'}W[y+\mu(A^{-1}-S^{-1})e]}
#'
#' where \eqn{W=(A^{-2}+S^{-2})^{-1}}. This is a form of the primal log
#' barrier algorithm described above. Setting \eqn{\mu=0} in each step
#' yields an affine scaling variant of the algorithm. The basic linear
#' algebra of each iteration is essentially unchanged, only the form
#' of the diagonal weighting matrix \eqn{W} has chagned.
#'
#' @references Portnoy S, Koenker R. The Gaussian hare and the Laplacian tortoise:
#' computability of squared-error versus absolute-error estimators.
#' \emph{Statistical Science}, 1997, 12(4): 279-300.
#' @export
#' @author Wenjing Wang \email{wenjingwangr@gmail.com}
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(quantreg)
#' data(ais)
#' tau <- c(0.1, 0.5, 0.9)
#' object <-rq(BMI ~ LBM + Ht, tau = tau, data = ais, method = 'fn')
#' frame_fn <- frame_fn_path(object, tau)
#' #plot the path
#' frame_fn1 <- frame_fn[[1]]
#' ggplot(frame_fn1, aes(x = value, y = objective)) +
#'    geom_point() +
#'    geom_path() +
#'    facet_wrap(~ variable, scale = 'free')
#'}
#'@useDynLib quokar

frame_fn_path <- function(object, tau){
    y <- matrix(object$model[, 1], ncol = 1)
    colnames(y) <- 'y'
    if(attr(object$terms, "intercept") == 1){
        x <- cbind(1, object$model[, -1])
    }else{
        x <- object$model[, -1]
    }
    z <- list()
  path <- list()
  path_m <- list()
  beta <- 0.99995
  eps <- 1e-06
  n <- length(y)
  x <- as.matrix(x)
  p <- ncol(x)
  if (n != nrow(x))
    stop("x and y don't match n")
  d <- rep(1, n)
  u <- rep(1, n)
  wn <- rep(0, 10 * n)
  1:length(tau) %>%
    map(function(i){
      rhs <- (1 - tau[i]) * apply(x, 2, sum)
      wn[1:n] <- (1 - tau[i])
      z[[i]] <- .Fortran("rqfnb",
                         n = as.integer(n),
                         p = as.integer(p),
                         a = as.double(t(as.matrix(x))),
                         c = as.double(-y),
                         rhs = as.double(rhs),
                         d = as.double(d),
                         u = as.double(u),
                         beta = as.double(beta),
                         eps = as.double(eps),
                         wn = as.double(wn),
                         wp = double((p + 3) * p),
                         it.count = integer(3),
                         info = integer(1),
                         it.routine = double(50*p),
                         PACKAGE = "quokar")
      iter <- z[[i]]$it.count[1] + 1
      coef_path <- -matrix(z[[i]]$it.routine, ncol = p,
                           byrow = TRUE)[1:iter, ]
      nn <- nrow(coef_path)
      objective <- rep(0, nn)
      for (j in 1:nn){
        resid <- data.frame(resid = y - x %*%
                              matrix(coef_path[j, ], ncol = 1))
        left_half <- abs(resid[resid <= 0])
        right_half <- abs(resid[resid > 0])
        m1 <- (1 - tau[i]) * left_half
        m2 <- tau[i] * right_half
        objective[j] <- sum(m1) + sum(m2)
      }
           if(attr(object$terms, 'intercept') == 1){
          coef_path <- -matrix(z[[i]]$it.routine, ncol = p,
                               byrow = TRUE)[1:iter, -1]
      }else{
          coef_path <- -matrix(z[[i]]$it.routine, ncol = p,
                               byrow = TRUE)[1:iter, ]
      }
      colnames(coef_path) <- paste("x", 1:ncol(coef_path), sep="")
      path[[i]] <- data.frame(objective, coef_path)
      path_m[[i]] <- gather(path[[i]], variable,
                                   value, -objective)
      tau_flag <- paste('tau=', tau[i], sep = "")
      path_m[[i]] <- cbind(tau_flag, path_m[[i]])
    })
}
