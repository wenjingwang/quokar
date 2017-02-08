globalVariables(c("variable", "value"))
#' Function to create data frame to plot the fitting path of quantile
#' regression using interior point method
#' @param y dependent variable in quantile regression
#' @param x design matrix in quantile regression
#' @param tau quantiles
#' @details This function used to illustrate the fitting process of
#' quantile regression using interior point method.
#' Koenker and Bassett(1978) introduced asymmetric weight on positive
#' and negative residuals, and solves the slightly modified l1-problem.
#'
#' $$min_{b \in R^{p}}\sum_{i=1}^{n}\rho_{\tau}(y_i-x_{i}^{'}b)$$
#'
#' where $\rho_{\tau}(r)=r[\tau-I(r<0)]$ for $\tau \in (0,1)$. This
#' yields the modified linear program
#'
#' $$min(\tau e^{'}u+(1-\tau)e^{'}v|y=Xb+u-v,(u,v) \in R_{+}^{2n})$$
#'
#' Adding slack variables, $s$, satisfying the constrains $a+s=e$, we
#' obtain the barrier function
#'
#' $$B(a, s, u) = y^{'}a+\mu \sum_{i=1}^{n}(loga_{i}+logs_{i})$$
#'
#' which should be maximized subject to the constrains
#' $X^{'}a=(1-\tau)X^{'}e$ and $a+s=e$. The Newton step $\delta_{a}$
#' solving
#'
#' $$max{y^{'}\delta_{a}+\mu \delta^{'}_{a}(A^{-1}-S^{-1})e-
#' \frac{1}{2}\mu \delta_{a}^{'}(A^{-2}+S^{-2})\delta_{a}}$$
#'
#' subject to $X{'}\delta_{a}=0$, satisfies
#'
#' $$y+\mu(A^{-1}-S^{-1})e-\mu(A^{-2}+S^{-2})\delta_{a}=Xb$$
#'
#' for some $b\in R^{p}$, and \delta_{a}$ such that $X^{'}\delta_{a}=0$.
#' Using the constraint, we can solve explicitly for the vector $b$,
#'
#' $$b=(X^{'}WX)^{-1}X^{'}W[y+\mu (A^{-1}-S^{-1})e]$$
#'
#' where $W=(A^{-2}+S^{-2})^{-1}$. This is a form of the primal log
#' barrier algorithm described above. Setting $\mu=0$ in each step
#' yields an affine scaling variant of the algorithm. The basic linear
#' algebra of each iteration is essentially unchanged, only the form
#' of the diagonal weighting matrix $W$ has chagned.
#'
#' @importFrom tidyr gather, ggplot2 ggplot
#' @references Portnoy S, Koenker R. The Gaussian hare and the Laplacian tortoise:
#' computability of squared-error versus absolute-error estimators[J].
#' Statistical Science, 1997, 12(4): 279-300.
#' @export
#'
frame_fn <- function(y, x, tau){
  z <- list()
  path <- list()
  path_m <- list()
  1:length(tau) %>%
    map(function(i){
      z[[i]] <- .Fortran(.tau[i]...)
      coef_path <- z[[i]]$routine
      coef_path <- as.matrix(coef_path)
      colnames(coef_path) <- paste('x', 1:ncol(coef_path), sep = '')
      n <- nrow(coef_path)
      objective <- rep(0, n)
      for (j in 1:n){
        resid <- data.frame(resid = y - x %*%
                              matrix(coef_path[j, ], ncol = 1))
        left_half <- resid[resid <= 0]
        right_half <- resid[resid > 0]
        m1 <- (1 - tau[i]) * left_half
        m2 <- tau[i] * right_half
        objective[j] <- sum(m1) + sum(m2)
      }
      path[[i]] <- data.frame(objective, coef_path)
      path_m[[i]] <- path[[i]] %>%
        gather(variable, value, -objective)
      tau_flag <- paste('tau=', tau[i], sep = "")
      path_m[[i]] <- cbind(tau_flag, path_m[[i]])
    })


  z <- list(z1 = list(routine = matrix(1:20, ncol = 2)),
            z2 = list(routine = matrix(2:21, ncol = 2)),
            z3 = list(routine = matrix(3:22, ncol = 2)))
  y <- matrix(3:22, ncol = 1)
  x <- cbind(1, 4:23)
  nrow(x)
  tau <- c(0.1, 0.5, 0.9)
