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

ais_bayes <- bayesQR(y ~ x - 1 , quantile = 0.5, alasso = FALSE, ndraw = 10000, prior = NULL)
ais_bayes_beta <- summary(ais_bayes)[[1]]$betadraw[, 1]
ais_bayes_sigma <- summary(ais_bayes)[[1]]$sigmadraw[, 1]

tau <- 0.5
M = 1000
n <- length(y)
beta <- ais_bayes_beta
sigma <- ais_bayes_sigma

KLD <- bayesKLD(y, x, beta, sigma, tau, M)

case <- 1:n
ais_bayes_data <- data.frame(case = case, KLD = KLD)

ggplot(ais_bayes_data, aes(x = case, y = KLD)) +
  ylab("Kullback-Leibler Divergence") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(ais_bayes_data, KLD > 0.35),
            aes(x = case, y = KLD-0.03, label = case))
