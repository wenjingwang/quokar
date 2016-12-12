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
x <- cbind(ais$LBM, sexInd)

tau <- 0.5
n <- length(y)
M = 1000
beta <- ais_bayes_beta
sigma <- ais_bayes_sigma

po <- bayesProb(y, x, M, method = "bayes.prob")

case <- 1:202
ais_bayes_data <- data.frame(case = case, po = po)

ggplot(ais_bayes_data, aes(x = case, y = po)) +
  ylab("posterior probability of being an outlier") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(ais_bayes_data, po > 0.004),
            aes(x = case, y = po-0.0002, label = case))
