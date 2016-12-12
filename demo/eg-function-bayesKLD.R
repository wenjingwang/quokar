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
x <- cbind(ais$LBM, sexInd)

tau <- 0.5
M = 1000
n <- length(y)

KLD <- qrod_bayes(y, x, tau, M, method = "bayes.kl")

case <- 1:n
ais_bayes_data <- data.frame(case = case, KLD = KLD)

ggplot(ais_bayes_data, aes(x = case, y = KLD)) +
  ylab("Kullback-Leibler Divergence") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(ais_bayes_data, KLD > 0.35),
            aes(x = case, y = KLD-0.03, label = case))
