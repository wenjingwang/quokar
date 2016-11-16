ALDqr_GCD <- function(y, x = NLLL, tau = NULL,
                      error = 1e-06, iter = 2000)
{
  p <- ncol(x)
  n <- length(y)
  theta <- EM.qr(y, x, tau, error, iter)$theta
  beta_qr <- theta[1:p, ]
  sigma_qr <- theta[p+1]
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
  gamma2 <- (2 + thep^2/taup2)/sigma_qr
  muc <- y - x %*% beta_qr
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                  gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)

  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                  gamma2), 0.5)) * (sqrt(delta2 / gamma2))
  E1 <- matrix(0, nrow = p, ncol = n)
  for(i in 1:n){
    suma2 <- x[-i,] * c(vchpN[-i] * (y[-i] - x[-i,] %*% beta_qr) - thep)
    E1[,i] <- apply(suma2, 2, sum)/(taup2)
  }
  E2 <- 1: n %>%
    map(function(i) {
      muc_i <- y[-i] - x[-i, ]%*%beta_qr
      sum(3*sigma_qr - (vchpN[-i] * muc_i^2 -
                          2 * muc_i * thep + vchp1[-i] *(thep^2 + 2 * taup2))/taup2)
    })
  E2 <- simplify2array(E2)
  Q1_beta <- E1/sigma_qr
  Q1_sigma <- -E2/(2*sigma_qr^2)
  xM <- c(sqrt(vchpN)) * x
  suma1 <- t(xM) %*% (xM)
  Q2_beta <- -(suma1)/(sigma_qr * taup2)
  Q2_sigma <- 3/(2*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +
                                        vchp1*(thep^2 + 2*taup2)))/(sigma_qr^3*taup2)
  GCD_beta <- 1:n %>%
    map(function(i){
      c(Q1_beta[,i]) %*% solve(-Q2_beta) %*% matrix(Q1_beta[,i], ncol = 1)
    })
  GCD_beta <- simplify2array(GCD_beta)
  GCD_sigma <- 1:n %>%
    map(function(i) {
      Q1_sigma[i]*solve(-Q2_sigma)*Q1_sigma[i]
    })
  GCD_sigma <- simplify2array(GCD_sigma)
  GCD <- GCD_beta + GCD_sigma
  return(GCD)
}
