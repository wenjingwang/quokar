#Kullback-Leibler Divergence

library(bayesQR)
library(MCMCpack)
library(coda)
library(MASS)
library(ALDqr)
library(GIGrvg)
library(ggplot2)
library(caTools)


data(ais)
y <- ais$BMI
sexInd <-(ais$Sex == 1) + 0
x <- cbind(1, ais$LBM, sexInd)

ais_bayes <- bayesQR(y ~ x -1 , quantile = 0.5, alasso = FALSE, ndraw = 10000, prior = NULL)
ais_bayes_beta <- summary(ais_bayes)[[1]]$betadraw[, 1]
ais_bayes_sigma <- summary(ais_bayes)[[1]]$sigmadraw[, 1]


tau <- 0.5
beta <- ais_bayes_beta
sigma <- ais_bayes_sigma
taup2 <- (2/(tau * (1 - tau)))
theta <- (1 - 2 * tau) / (tau * (1 - tau))
n <- length(y)

delta <- sqrt((y - x %*% beta)^2/(sigma * taup2))
gamma <- sqrt(theta^2/(sigma*taup2) + 2/sigma)

M <- 1000
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

case <- 1:202
ais_bayes_data <- data.frame(case = case, KLD = KLD)

ggplot(ais_bayes_data, aes(x = case, y = KLD)) +
  ylab("Kullback-Leibler Divergence") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(ais_bayes_data, KLD > 0.35),
            aes(x = case, y = KLD-0.03, label = case))
