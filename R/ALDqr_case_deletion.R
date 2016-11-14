ALDqr_case_deletion <- function(y, x = NLLL, theta, sigma) {

  p <- ncol(x)

  n <- nrow(x)

  taup2 <- (2/(tau * (1 - tau)))

  thep <- (1 - 2 * tau) / (tau * (1 - tau))

  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)

  gamma2 <- (2 + thep^2/taup2)/sigma_qr

  muc <- y - x %*% beta_qr

  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
          gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)

  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
          gamma2), 0.5)) * (sqrt(delta2 / gamma2))

  E1 <- 1: length(y) %>%
    map(function(i){
      E1_i <- apply((vchpN[-i]%*%(y[-i]-x[-i,]%*% beta_qr))[,]*x[-i,]- x[-i,]*thep, 2,                           sum)/taup2
      return(E1_i)
    })

  E1 <- simplify2array(E1)

  E2 <- 1: length(y) %>%
    map(function(i) {
      E2_i <- sum(matrix(3*sigma_qr*c(rep(1,n-1)), ncol = 1)- vchpN[-i,]%*%
                    (t(y[-i]-x[-i,]%*%beta_qr)%*%(y[-i]-x[-i,] %*%
                  beta_qr))- 2*(y[-i]-x[-i,]%*% beta_qr)*thep +
                    matrix(vchp1[-i,]*thep^2/4, ncol = 1))
      return(E2_i)
    })

  E2 <- simplify2array(E2)

  Q2_sigma <- 3/(4*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +vchp1*(thep^2 + 2*taup2)))/(2*sigma_qr^3*taup2)

  xM <- c(sqrt(vchpN)) * x

  suma1 <- t(xM) %*% (xM)

  beta_i <- beta_qr + taup2*solve(suma1)%*%E1[[1]]

  sigma_i <- sigma_qr^2 - 1/(Q2_sigma)*E2[[1]]/(2*sigma_qr^2)
}
