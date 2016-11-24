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
n <- length(y)
M = 1000
beta <- ais_bayes_beta
sigma <- ais_bayes_sigma

po <- bayesProb(y, x, beta, sigma, tau, M)

case <- 1:202
ais_bayes_data <- data.frame(case = case, po = po)

ggplot(ais_bayes_data, aes(x = case, y = po)) +
  ylab("posterior probability of being an outlier") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(ais_bayes_data, po > 0.004),
            aes(x = case, y = po-0.0002, label = case))
