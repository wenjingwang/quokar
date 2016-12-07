#'Outlier Dignostic for quantile regression
#'@para y dependent variable in quantile regression
#'
#'@para x indepdent variables in quantile regression.
#'Note that: x is the independent variable matrix which including
#'the intercept. That means, if the dimension of independent
#'variables is p and the sample size is n, x is a (n \times (p+1))
#'matrix with the first column is 1.
#'
#'@para error The EM algorithm accuracy of error used in MLE estimation
#'
#'@para iter the iteration frequancy for EM algorithm used in MLE estimation
#'
#'@M the iteration frequancy for MCMC used in Baysian Estimation
#'
#'@description
#'This group of function is used to compute diagnositcs for a
#'quantile regression model based on the building blocks returned
#'by.
#'
#'@details
#'The primary function is which returns either
#'a list or data frame of influence mesures depending on whether
#'or if only one aspect of the model is selected
#'then a list with Generalized Cook's distance,Q-function Distance,
#'Baysian Method or ES method.
#'
#'This function is used to detect outlier in quantile regression
#'level. It use code to fit the models for each quantile.
#'
#'The method cook.distance, qfunction.distance,
#'baysian.probability, baysian.distribution and ES
#'can be used for direct computation of the corresponding diagnostic
#'quantities from an object of class case_delete
#'
#'@author Wenjing Wang \email{wenjingwang1990@gmail.com}
#'@keywords quantile regression outlier diagnostics
#'@references Benites L E, Lachos V H, Vilca F E.(2015)``Case-Deletion
#' Diagnostics for Quantile Regression Using the Asymmetric Laplace
#' Distribution,\emph{arXiv preprint arXiv:1509.05099}.
#'
#'Hawkins D M, Bradu D, Kass G V.(1984)``Location of several outliers in
#'multiple-regression data using elemental sets. \emph{Technometrics},
#'\bold{26(3)}, 197-208.
#'
#'Koenker R, Bassett Jr G.(1978)`` Regression quantiles,
#'\emph{Econometrica},\bold{1}, 33-50.
#'
#'Santos B, Bolfarine H.(2016)``On Baysian quantile regression and
#'outliers,\emph{arXiv:1601.07344}
#'
#'Kozumi H, Kobayashi G.(2011)``Gibbs sampling methods for Bayesian
#'quantile regression,\emph{Journal of statistical computation and
#'simulation}, \bold{81(11)}, 1565-1578.
#'
#'@examples
#'data(ais, package = `ALDqr`)
#'
#'
#'
#'

qrod <- function(y, x, tau, error = NULL, iter = NULL,
                 method = c("cook.distance", "qfunciton", "bayes.probability",
                            "bayes.KLD"), M = NULL){
  method <- match.arg(method)
    if(method == "cook.distance"){
      ALDqr_GCD <- function(y, x, tau, error, iter){
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
      distance <- ALDqr_GCD(y, x, tau, error, iter)
    }else if(method == "qfunction"){
      ALDqr_QD <- function(y, x, tau, error, iter){
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
      distance <- ALDqr_QD(y, x = NULL, tau = NULL, error, inter)
    }else if(method == "bayes.probability"){
      bayesKLD <- function(y, x, beta, sigma, tau, M){
        taup2 <- (2/(tau * (1 - tau)))
        theta <- (1 - 2 * tau) / (tau * (1 - tau))
        n <- length(y)
        delta <- sqrt((y - x %*% beta)^2/(sigma * taup2))
        gamma <- sqrt(theta^2/(sigma*taup2) + 2/sigma)
        nu_dens <- matrix(0, nrow = M, ncol = n)
        for(i in 1:n) {
          nu_dens[,i] <- rgig(M, 1/2, delta[i], gamma)
        }
        KLS <- matrix(0, nrow = n, ncol = n-1)
        for(i in 1:n) {
          nu1 <- nu_dens[, -i]
          hi <- density(nu_dens[, i], kernel = "gaussian")$bw
          upper_x <- max(nu_dens[, i])
          lower_x <- min(nu_dens[, i])
          for(j in 1:(n-1)) {
            hj <- density(nu1[, j], kernel = "gaussian")$bw
            func <- function(xx) {
              log((1/(M*sqrt(2*pi))*sum(exp(-1/2*((xx-nu_dens[,i])/hi)^2)))/
                    (1/(M*sqrt(2*pi))*sum(exp(-1/2*((xx-nu1[,j])/hj)^2))))*
                (1/(M*sqrt(2*pi))*sum(exp(-1/2*((xx-nu_dens[,i])/hi)^2)))
            }
            KLS[i,j] <- (upper_x - lower_x)*(func(upper_x)+func(lower_x))/2
          }
        }

        KLD <- apply(KLS, 1, mean)
        for(i in 1:n) {
          if (KLD[i] == Inf) {KLD[i] = 1}
        }
        return(KLD)
      }
      distance <- bayesProb(y, x, beta, sigma, tau, M)
    }else if(method == "bayes.KLD"){
      bayesProb <- function(y, x, beta, sigma, tau, M){
        taup2 <- (2/(tau * (1 - tau)))
        theta <- (1 - 2 * tau) / (tau * (1 - tau))
        delta <- sqrt((y - x %*% beta)^2/(sigma * taup2))
        gamma <- sqrt(theta^2/(sigma*taup2) + 2/sigma)
        nu_dens <- matrix(0, nrow = M, ncol = n)
        for(i in 1:n) {
          nu_dens[,i] <- rgig(M, 1/2, delta[i], gamma)
        }
        A <- matrix(0, nrow = n, ncol = n-1)
        for(i in 1:n) {
          probs <- 1:M/M
          nu1 <- nu_dens[, -i]
          for(j in 1:(n-1)) {
            A[i,j] <- 1/M*sum(quantile(nu_dens[,i], probs = probs) > max(nu1[, j]))
          }
        }
        prob <- apply(A, 1, mean)
        return(prob)
      }
      distance <- bayesKLD(y, x, beta, sigma, tau, M)
    }else if(method %in% c("cook.distance", "qfunction", "bayes.probability",
                            "bayes.KLD") == FALSE)
      warning("Method should be one of 'cook.distance', 'qfunction',
              'bayes.probability', 'bayes.KLD'")
  returen(distance)
}

