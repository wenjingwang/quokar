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
tau = 0.2
M = 1000

po <- qrod_bayes(y, x, tau, M, method = "bayes.prob")

ggplot(po, aes(x = case, y = result)) +
  ylab("posterior probability of being an outlier") +
  xlab("case") +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(po, result > 0.004),
            aes(x = case, y = result - 0.0002, label = case))
