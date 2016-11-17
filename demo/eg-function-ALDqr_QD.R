library(ALDqr)
library(ggplot2)
library(magrittr)
library(purrr)
library(gridExtra)

data(ais)
y <- ais$BMI
sexInd <-(ais$Sex == 1) + 0
x <- cbind(1, ais$LBM, sexInd)

tau <- c(0.1, 0.25, 0.5, 0.9)
n <- length(y)
case <- 1: n

QD <- ALDqr_QD(y, x, tau[1], error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p1 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 1),
            aes(case, QD+0.1, label = case))

QD <- ALDqr_QD(y, x, tau[2], error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p2 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 1),
            aes(case, QD+0.1, label = case))

QD <- ALDqr_QD(y, x, tau[3], error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p3 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 0.5),
            aes(case, QD+0.1, label = case))

QD <- ALDqr_QD(y, x, tau[4], error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p4 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 1),
            aes(case, QD+0.1, label = case))

grid.arrange(p1, p2, p3, p4, ncol = 2)
