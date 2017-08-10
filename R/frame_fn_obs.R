#' @title Visualization of quantile regression model
#' fitting: interior point algorithm
#' @param object quantile regression model using interior point method for estimating
#' @param tau quantile
#' @return Weighted observations in quantile regression fitting using interior point
#' algorithm
#' @details This function used to illustrate data used in fitting process of
#' quantile regression based on interior point method.
#' Koenker and Bassett(1978) introduced asymmetric weight on positive
#' and negative residuals, and solves the slightly modified l1-problem.
#'
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
#'
#' @export
#' @useDynLib quokar
#' @examples
#' library(ggplot2)
#' library(quantreg)
#' library(tidyr)
#' library(dplyr)
#' library(gridExtra)
#' data(ais)
#' tau <- c(0.1, 0.5, 0.9)
#' object <- rq(BMI ~ LBM + Ht, data = ais, tau = tau, method = 'fn')
#' fn <- frame_fn_obs(object, tau)
#' ##For tau = 0.1, plot the observations used in quantile regression
#' ##fitting based on interior point method
#' fn1 <- fn[ ,1]
#' case <- 1:length(fn1)
#' fn1 <- cbind(case, fn1)
#' m <- data.frame(y = ais$BMI, x1 = ais$LBM, x2 = ais$Ht, fn1)
#' p <- length(attr(object$coefficients, "dimnames")[[1]])
#' m_f <- m %>% gather(variable, value, -case, -fn1, -y)
#' mf_a <- m_f %>%
#'  group_by(variable) %>%
#'  arrange(variable, desc(fn1)) %>%
#'  filter(row_number() %in% 1:p )
#' p1 <- ggplot(m_f, aes(x = value, y = y)) +
#'  geom_point(alpha = 0.1) +
#'  geom_point(data = mf_a, size = 3) +
#'  facet_wrap(~variable, scale = "free_x")
#' ## For tau = 0.5, plot the observations used in quantile regression
#' ##fitting based on interior point method
#' fn2 <- fn[,2]
#' case <- 1: length(fn2)
#' fn2 <- cbind(case, fn2)
#' m <- data.frame(y = ais$BMI, x1 = ais$LBM, x2 = ais$Ht, fn2)
#' p <- length(attr(object$coefficients, "dimnames")[[1]])
#' m_f <- m %>% gather(variable, value, -case, -fn2, -y)
#' mf_a <- m_f %>%
#'    group_by(variable) %>%
#'    arrange(variable, desc(fn2)) %>%
#'    filter(row_number() %in% 1:p )
#' p2 <- ggplot(m_f, aes(x = value, y = y)) +
#'    geom_point(alpha = 0.1) +
#'    geom_point(data = mf_a, size = 3) +
#'    facet_wrap(~variable, scale = "free_x")
#' ## For tau = 0.9
#' fn3 <- fn[,3]
#' case <- 1: length(fn3)
#' fn3 <- cbind(case, fn3)
#' m <- data.frame(y = ais$BMI, x1 = ais$LBM, x2 = ais$Ht, fn3)
#' p <- length(attr(object$coefficients, "dimnames")[[1]])
#' m_f <- m %>% gather(variable, value, -case, -fn3, -y)
#' mf_a <- m_f %>%
#'   group_by(variable) %>%
#'   arrange(variable, desc(fn3)) %>%
#'   filter(row_number() %in% 1:p )
#' p3 <- ggplot(m_f, aes(x = value, y = y)) +
#'   geom_point(alpha = 0.1) +
#'   geom_point(data = mf_a, size = 3) +
#'   facet_wrap(~variable, scale = "free_x")
#' grid.arrange(p1, p2, p3, ncol = 1)

frame_fn_obs <- function(object, tau){
    y <- matrix(object$model[,1], ncol = 1)
    if(attr(object$terms, "intercept") == 1){
    x <- cbind(1, object$model[, -1])
    }else{
        x <- object$model[, -1]
    }
    x <- as.matrix(x)
    colnames(y) <-'y'
    z <- list()
    w_obs <- list()
    beta <- 0.99995
    eps <- 1e-06
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x)) stop("x and y don't match n")
    d <- rep(1, n)
    u <- rep(1, n)
    wn <- rep(0, 10 * n)
    ntau <- length(tau)
    W <- matrix(0, nrow = n, ncol = ntau)
    for(i in 1:ntau){
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
            W[, i] <- z[[i]]$d
            W[, i] <- z[[i]]$d/sum(z[[i]]$d)
    }
            #x <- as.matrix(object$model[, -1])
            #wx <- diag(sqrt(W)) %*% x
            #wy <- diag(sqrt(W)) %*% y
            #colnames(wx) <- paste('x', 1:ncol(x), sep='')
            #colnames(wy) <- 'y'
            #tau_flag <- paste('tau', tau[i], sep = '')
            #w_obs[[i]] <- data.frame(wy, wx)
            #w_obs[[i]] <- cbind(tau_flag, w_obs[[i]])
    colnames(W) <- paste("tau", tau, sep ="")
    return(W)
}



