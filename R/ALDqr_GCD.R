ALDqr_GCD <- function(y, x = NLLL, tau = NULL,
               error = 1e-06, iter = 2000)
{

  theta <- EM.qr(y, x, tau, error, iter)$theta

  theta_i <- ALDqr_case_deletion(y, x, tau, error, iter)

  beta_i <- theta_i$beta_i

  sigma_i <- theta_i$sigma_i

  xM <- c(sqrt(vchpN)) * x

  suma1 <- t(xM) %*% (xM)

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

  Q2_beta <- -solve(suma1)/ (sigma_qr * taup2)

  Q2_sigma <- 3/(4*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +
              vchp1*(thep^2 + 2*taup2)))/(2*sigma_qr^3*taup2)

  GcD <- rep(0, length(y))

1:length(y) %>%
  map(function(i){
    GcD[i] <- t(beta_i[[i]])%*%Q2_beta%*%beta_i[[i]] + sigma_i[i]*Q2_sigma*sigma_i[i]
    })

  return(GCD)
}


