library(ALDqr)
library(ggplot2)
library(magrittr)
library(purrr)
library(gridExtra)


data(ais)
y <- ais$BMI
sexInd <-(ais$Sex == 1) + 0
x <- cbind(1, ais$LBM, sexInd)

tau <- c(0.1, 0.2, 0.5, 0.9)
n <- length(y)
case <- 1: n

plots <- list()
for(i in seq_along(tau)){
  QD <- ALDqr_QD(y, x, tau[i], error = 1e-06, iter = 2000)
  data_QD <- data.frame(case = case, QD = QD)

  plots[[i]] <- ggplot(data_QD, aes(x = case, y = QD)) +
    geom_point(colour = "darkblue") +
    geom_text(data = subset(data_QD, QD > 1),
              aes(case, QD+0.1, label = case))
}
grid.arrange(grobs=plots)

