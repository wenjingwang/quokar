#Kullback-Leibler Divergence
library(bayesQR)
library(MCMCpack)
library(coda)
library(MASS)
library(ALDqr)
library(GIGrvg)
library(ggplot2)
library(caTools)
library(quokar)
data(ais)
y <- ais$BMI
sexInd <-(ais$Sex == 1) + 0
x <- cbind(ais$LBM, sexInd)

tau <- 0.5
M = 1000
n <- length(y)

KLD <- qrod_bayes(y, x, tau, M, method = "bayes.kl")

ggplot(KLD, aes(x = case, y = result)) +
  ylab("Kullback-Leibler Divergence") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(KLD, result > 0.35),
            aes(x = case, y = result - 0.03, label = case))
