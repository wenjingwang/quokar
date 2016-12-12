data(ais)
y <- ais$BMI
sexInd <-(ais$Sex == 1) + 0
x <- cbind(1, ais$LBM, sexInd)

tau <- c(0.1, 0.25, 0.5, 0.9)
n <- length(y)
case <- 1: n

plots <- list()
for(i in seq_along(tau)){
  GCD <- qrod_mle(y, x, tau[i], error = 1e-06, iter = 2000, method = 'qfunction')
  data_GCD <- data.frame(case = case, GCD = GCD)
  critical_GCD <- 2*(tau[i] + 1)/n
  plots[[i]] <- ggplot(data_GCD, aes(x = case, y = GCD)) +
    geom_point(colour = "darkblue") +
    geom_hline(yintercept = critical_GCD, colour = "red") +
    geom_text(data = subset(data_GCD, GCD > critical_GCD),
              aes(case, GCD + 0.001, label = case))
}
grid.arrange(grobs=plots)
