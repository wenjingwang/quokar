library(ALDqr)
library(ggplot2)
library(magrittr)
library(purrr)
data(ais)
y <- ais$BMI
sexInd <-(ais$Sex == 1) + 0
x <- cbind(1, ais$LBM, sexInd)
#Calculate Generalized Cook Distance
tau <- c(0.1, 0.25, 0.5, 0.9)
for (i in 1:4) {
  i=0.9
  tau <- tau[i]
  GCD <- ALDqr_GCD(y, x, tau, error = 1e-06, iter = 2000)
  n <- length(y)
  case <- 1: n
  data_GCD <- data.frame(case = case, GCD = GCD)
  critical_GCD <- 2*(tau + 1)/n
  ggplot(data_GCD, aes(x = case, y = GCD)) +
          geom_point(colour = "darkblue") +
          geom_hline(yintercept = critical_GCD, colour = "red") +
          geom_text(data = subset(data_GCD, GCD > critical_GCD),
            aes(case, GCD + 0.001, label = case))
}
