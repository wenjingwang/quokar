library(ggplot2)
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
  QD <- qrod_mle(y, x, tau[i], error = 1e-06, iter = 2000, method = "cook.distance")
  plots[[i]] <- ggplot(QD, aes(x = case, y = distance)) +
    geom_point(colour = "darkblue") +
    geom_text(data = subset(QD, distance > 0.001),
              aes(case, distance, label = case))
}
grid.arrange(grobs=plots)

