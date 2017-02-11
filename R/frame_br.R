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
#' @importFrom purrr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr inner_join
#' @importFrom stats qt qnorm lm coefficients
#' @importFrom quantreg bandwidth.rq rq.fit.br
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
#'
frame_br <- function(object, tau){
  wh <- function(object, tau){
    x <- object$x
    y <- object$y
    tol <- .Machine$double.eps^(2/3)
    eps <- tol
    big <- .Machine$double.xmax
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    ny <- NCOL(y)
    nsol <- 2
    ndsol <- 2
    if (qr(x)$rank < p)
      stop("Singular design matrix")
    if (tau < 0 || tau > 1) {
      nsol <- 3 * n
      ndsol <- 3 * n
      lci1 <- FALSE
      qn <- rep(0, p)
      cutoff <- 0
      tau <- -1
    }
    else {
      if (p == 1)
        ci <- FALSE
      if (ci) {
        lci1 <- TRUE
        if (tcrit)
          cutoff <- qt(1 - alpha/2, n - p)
        else cutoff <- qnorm(1 - alpha/2)
        if (!iid) {
          h <- bandwidth.rq(tau, n, hs = TRUE)
          bhi <- rq.fit.br(x, y, tau + h, ci = FALSE)
          bhi <- coefficients(bhi)
          blo <- rq.fit.br(x, y, tau - h, ci = FALSE)
          blo <- coefficients(blo)
          dyhat <- x %*% (bhi - blo)
          if (any(dyhat <= 0)) {
            pfis <- (100 * sum(dyhat <= 0))/n
            warning(paste(pfis, "percent fis <=0"))
          }
          f <- pmax(eps, (2 * h)/(dyhat - eps))
          qn <- rep(0, p)
          for (j in 1:p) {
            qnj <- lm(x[, j] ~ x[, -j] - 1, weights = f)$resid
            qn[j] <- sum(qnj * qnj)
          }
        }
        else qn <- 1/diag(solve(crossprod(x)))
      }
      else {
        lci1 <- FALSE
        qn <- rep(0, p)
        cutoff <- 0
      }
    }
    z <- .Fortran("rqbr",
                  n = as.integer(n),
                  p = as.integer(p),
                  n5 = as.integer( n + 5),
                  p3 = as.integer(p + 3),
                  p4 = as.integer(p + 4),
                  a = as.double(x),
                  b = as.double(y),
                  t = as.double(tau),
                  tole = as.double(tol),
                  flag = as.integer(1),
                  coef = double(p),
                  resid = double(n),
                  s =  integer(n),
                  wa = double((n + 5) * (p + 4)),
                  wb = double(n),
                  nsol = as.integer(nsol),
                  ndsol = as.integer(ndsol),
                  sol = double((p + 3) * nsol),
                  dsol = double(n * ndsol),
                  lsol = as.integer(0),
                  wh = integer(p * nsol),
                  qn = as.double(qn),
                  cutoff = as.double(cutoff),
                  ci = double(4 * p),
                  tnmat = double(4 *  p),
                  big =  as.double(big),
                  lci1 = as.logical(lci1),
                  PACKAGE = "quokar")
    return(z$wh[1:NCOL(x)])
  }
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
