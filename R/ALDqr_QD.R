ALDqr_QD <- function(y, x = NLLL, tau = NULL,
                     error = 1e-06, iter = 2000)
{
  p <- ncol(x)
  n <- length(y)
  theta_all <- EM.qr(y, x, tau, error, iter)$theta
  beta_all <- theta_all[1:p, ]
  sigma_all <- theta_all[p+1]
  beta_i <- ALDqr_case_deletion(y, x, tau, error, iter)$beta_i
  sigma_i <-  ALDqr_case_deletion(y, x, tau, error, iter)$sigma_i
  Q_function <- function(beta_qr, sigma_qr, tau){
    n <- length(y)
    taup2 <- (2/(tau * (1 - tau)))
    thep <- (1 - 2 * tau) / (tau * (1 - tau))
    delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
    gamma2 <- (2 + thep^2/taup2)/sigma_qr
    muc <- y - x %*% beta_qr
    vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                    gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)
    vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                    gamma2), 0.5)) * (sqrt(delta2 / gamma2))
    Q <- (-3*log(sigma_qr)/2)*n - sum(vchpN * muc^2 - 2 * muc * thep +
                vchp1 *(thep^2 + 2 * taup2))/(2 * sigma_qr * taup2)
    return(Q)
  }
  Q_all <- Q_function(beta_all, sigma_all, tau)
  Q_i <- rep(0, n)
  for(i in 1:n){
    Q_i[i] <- Q_function(beta_i[,i], sigma_i[i], tau)
  }
  QD <- 1:n %>%
    map(function(i) 2*(Q_all - Q_i[i]))
  QD <- simplify2array(QD)
  return(QD)
}
