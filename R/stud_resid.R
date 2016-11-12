#' studentised residuals
#'
#' Function for computing studentised residuals
#' from a quantile regression fit
#' @param object an object of class quantile regression model
#' @details The \code{stud_resid} function provides a studentized
#' residual on each quantile.
#' @author Wenjing Wang \email{wenjingwang1990@gmail.com}
#' @references
#' Ranganai, E(2016). On studentized residuals in the quantile
#' regression framework. University of South Africa.
#' \emh{SpringerPlus}, \bold{5}(1), 1231.


stud_resid <- function(object) {

  xj <- cbind(1, na.exclude(es(object)))

  xi <- object$x[unclass(attr(na.exclude(es(object)), "na.action")), ]

  lm_coef <- lm(object$y[unclass(attr(na.exclude(es(object)), "na.action"))] ~ xi[, 2:ncol(xi)])$coef

  n = length(lm_coef)

  error <- object$y - (object$x %*% matrix(lm_coef, nrow = n))

  hii <- diag(object$x %*% (solve(t(xj) %*% xj)) %*% t(object$x))

  epson <- error / sqrt(1 + hii)
  epson2 <- epson^2

  press_i <- simplify2array( 1:length(object$y) %>%
    map(function(i) {sum(epson2) -  epson2[i]}))

  sepr <- epson / sqrt(press_i / (length(object$y) - 2*ncol(object$x) - 1))

  return(sepr)
}


stud_resid(model1)
