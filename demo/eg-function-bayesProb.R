library(bayesQR)
library(MCMCpack)
library(coda)
library(MASS)
library(ALDqr)
library(GIGrvg)
library(ggplot2)

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

A <- matrix(0, nrow = n, ncol = n-1)
for(i in 1:n) {
  probs <- 1:M/M
  nu1 <- nu_dens[, -i]
  for(j in 1:(n-1)) {
    A[i,j] <- 1/M*sum(quantile(nu_dens[,i], probs = probs) > max(nu1[, j]))
  }
}
po <- apply(A, 1, mean)
po
case <- 1:202
ais_bayes_data <- data.frame(case = case, po = po)

ggplot(ais_bayes_data, aes(x = case, y = po)) +
  ylab("posterior probability of being an outlier") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(ais_bayes_data, po > 0.004),
            aes(x = case, y = po-0.0002, label = case))
