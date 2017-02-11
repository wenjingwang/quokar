#globalVariables(c("obs", "r", "d"))
#' @title Density function plot for quantile regression fitting using
#' Asymmetric Laplace Distribution
#'
#' @param y vector, dependent variable of quantile regression
#' @param x matrix, matrix consisted independent variables of quantie
#' regression
#' @param tau sigle number or vector, quantiles
#' @param smooth sigular, default is 100, the larger the smoother of
#' density function
#' @param error the convergence maximum error
#' @param iter maximum iterations of the EM algorithm
#' @description density function used in quantile regression fitting
#' @export
#' @examples
#' library(ggplot2)
#' data(ais)
#' x <- matrix(ais$LBM, ncol = 1)
#' y <- ais$BMI
#' tau = c(0.1, 0.5, 0.9)
#' ald_data <- frame_ald(y, x, tau, smooth = 10, error = 1e-6,
#'                   iter = 2000)
#' ggplot(ald_data) +
#'    geom_line(aes(x = r, y = d, group = obs, colour = tau_flag)) +
#'    facet_wrap(~tau_flag, ncol = 1) +
#'    xlab('') +
#'    ylab('Asymmetric Laplace Distribution Density Function')
#'

frame_ald <- function(y, x, tau, smooth, error, iter){
  n <- length(y)
  ntau <- length(tau)
  y <- matrix(y, ncol = 1)
  colnames(y) <- 'y'
  x <- cbind(1, x)
  coef_qr <- list()
  beta_qr <- matrix(0, nrow = ncol(x), ncol = ntau)
  sigma_qr <- rep(0, ntau)
  xald <- matrix(0, nrow = n, ncol = ntau)
  for(i in 1:ntau){
    coef_qr[[i]] <- EM.qr(y, x, error, tau = tau[i],
                          iter, envelope = FALSE)
    beta_qr[,i] <- coef_qr[[i]]$theta[1:2, ]
    sigma_qr[i] <- coef_qr[[i]]$theta[3,]
    xald[,i] <- x %*% matrix(beta_qr[, i], ncol = 1)
  }
  rald <- matrix(0, nrow = smooth, ncol = n)
  dald <- matrix(0, nrow = smooth, ncol = n)
  rald_list <- matrix(0, nrow = smooth*ntau, ncol = n)
  dald_list <- matrix(0, nrow = smooth*ntau, ncol = n)
  tau_flag <- rep(0, smooth*ntau)
  m=1
  for(j in 1:ntau){
    p <- tau[j]
    for (i in 1:n){
      mu_i <- xald[i, j]
      sigma_i <- sigma_qr[j]
      rald[, i] <- rALD(smooth, mu_i, sigma_i, p)
      dald[, i] <- dALD(y = rald[,i], mu_i, sigma_i, p)
    }
    rald_list[m:(m+smooth-1), ] <- round(rald,2)
    dald_list[m:(m+smooth-1), ] <- round(dald,2)
    tau_flag[m:(m+smooth-1)] <- paste('tau=',
                                      rep(tau[j],smooth),
                                      sep = '')
    m = m + smooth
  }
  tau_flag <- data.frame(tau_flag)
  rald_list <- data.frame(rald_list)
  rald_list_m <- data.frame(tau_flag, rald_list)
  dald_list <- data.frame(dald_list)
  dald_list_m <- data.frame(tau_flag, dald_list)
  rald_list_g <- rald_list_m %>% gather(obs, r, -tau_flag)
  dald_list_g <- dald_list_m %>% gather(obs, d, -tau_flag)
  t_g <- data.frame(tau_flag = rald_list_g$tau_flag,
                    obs = rald_list_g$obs,
                    r = round(as.numeric(rald_list_g$r), 2),
                    d = round(as.numeric(dald_list_g$d), 2))
  t_g_a <- group_by(t_g, tau_flag) %>% arrange(r)
  return(t_g_a)
}
