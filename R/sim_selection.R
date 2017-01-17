library(quantreg)
library(quokar)
library(ggplot2)
library(gridExtra)
#simulation studies

#the errors are simulated from four different distrubtions

#a standard normal distribution
#a student-t distribution
#a heteroscedastic normal distribution
#a bimodal mixture distribution


##consider sample size
n = 50 # n = 50, 100, 150, 200, 300, 400, 500, 700, 800
beta1 = 1
beta2 = 1
beta3 = 1
beta <- cbind(beta1, beta2, beta3)

x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
x <- cbind(x1, x2, x3)

epsilon1 <- rnorm(n, 0, 1)
epsilon2 <- sgt::rsgt(n, 0, 2)
epsilon3 <- (1 + x2) * rnorm(n, 0, 2) #a heteroscedastic normal distribution
epsilon4 <- 0.6 * sgt::rsgt(n, -20, 2) + 0.4 * sgt::rsgt(n, 15, 2)

y1 <- x %*% t(beta) + epsilon1
y2 <- x %*% t(beta) + epsilon2
y3 <- x %*% t(beta) + epsilon3
y4 <- x %*% t(beta) + epsilon4
##data sets
data1 <- data.frame(y = y1, x1 = x1, x2 = x2, x3 = x3)
data2 <- data.frame(y = y2, x1 = x1, x2 = x2, x3 = x3)
data3 <- data.frame(y = y3, x1 = x1, x2 = x2, x3 = x3)
data4 <- data.frame(y = y4, x1 = x1, x2 = x2, x3 = x3)
##adding outliers
x_outlier <- matrix(runif(5*3, 0, 1), nrow = 5, ncol = 3)
y_outlier <- runif(5, 0, 1)
data_outlier <- data.frame(y = y_outlier, x1 = x_outlier[, 1],
                           x2 = x_outlier[, 2], x3 = x_outlier[, 3])
data1_outlier <- rbind(data1, data_outlier)
data2_outlier <- rbind(data2, data_outlier)
data3_outlier <- rbind(data3, data_outlier)
data4_outlier <- rbind(data4, data_outlier)

#######
data1_outlier <- cbind(case = 1:nrow(data1_outlier), data1_outlier)
ggplot(data1_outlier, aes(x = case, y = y)) +
  geom_point(colour = "darkblue")
###look at the simulation result##
x_model <- as.matrix(cbind(1, data1_outlier[, -1]))
method1 <- qrod_mle(y = data1_outlier[, 1], x = x_model, tau = 0.1, error = 1e-06,
                    iter = 1000, method = "cook.distance")
method2 <- qrod_mle(y = data1_outlier[, 1], x = x_model, tau = 0.1, error = 1e-06,
                    iter = 1000, method = "qfunction")
method3 <- qrod_bayes(y = data1_outlier[, 1], x = as.matrix(data1_outlier[, -1]),
                      tau = 0.1, M = 1000, method = "bayes.prob")
p1 <- ggplot(method1, aes(x = case, y = distance)) +
          geom_point(colour = "darkblue") +
          ggtitle("tau = 0.1")

p2 <- ggplot(method2, aes(x = case, y = distance)) +
          geom_point(colour = "darkblue") +
          ggtitle("tau = 0.1")

p3 <- ggplot(method3, aes(x = case, y = result)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.1")

method4 <- qrod_mle(y = data1_outlier[, 1], x = x_model, tau = 0.5, error = 1e-06,
                    iter = 1000, method = "cook.distance")
method5 <- qrod_mle(y = data1_outlier[, 1], x = x_model, tau = 0.5, error = 1e-06,
                    iter = 1000, method = "qfunction")
method6 <- qrod_bayes(y = data1_outlier[, 1], x = as.matrix(data1_outlier[, -1]),
                      tau = 0.5, M = 1000, method = "bayes.prob")

p4 <- ggplot(method4, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.5")

p5 <- ggplot(method5, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.5")

p6 <- ggplot(method6, aes(x = case, y = result)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.5")

method7 <- qrod_mle(y = data1_outlier[, 1], x = x_model, tau = 0.9, error = 1e-06,
                    iter = 1000, method = "cook.distance")
method8 <- qrod_mle(y = data1_outlier[, 1], x = x_model, tau = 0.9, error = 1e-06,
                    iter = 1000, method = "qfunction")
method9 <- qrod_bayes(y = data1_outlier[, 1], x = as.matrix(data1_outlier[, -1]),
                      tau = 0.9, M = 1000, method = "bayes.prob")
p7 <- ggplot(method7, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.9")

p8 <- ggplot(method8, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.9")

p9 <- ggplot(method9, aes(x = case, y = result)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.9")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)


##############################################################
data2_outlier <- cbind(case = 1:nrow(data2_outlier), data2_outlier)
ggplot(data2_outlier, aes(x = case, y = y)) +
  geom_point(colour = "darkblue")
###look at the simulation result##
x_model <- as.matrix(cbind(1, data2_outlier[, -1]))
method1 <- qrod_mle(y = data2_outlier[, 1], x = x_model, tau = 0.1, error = 1e-06,
                    iter = 1000, method = "cook.distance")
method2 <- qrod_mle(y = data2_outlier[, 1], x = x_model, tau = 0.1, error = 1e-06,
                    iter = 1000, method = "qfunction")
method3 <- qrod_bayes(y = data2_outlier[, 1], x = as.matrix(data2_outlier[, -1]),
                      tau = 0.1, M = 1000, method = "bayes.prob")
p1 <- ggplot(method1, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.1")

p2 <- ggplot(method2, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.1")

p3 <- ggplot(method3, aes(x = case, y = result)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.1")

method4 <- qrod_mle(y = data2_outlier[, 1], x = x_model, tau = 0.5, error = 1e-06,
                    iter = 1000, method = "cook.distance")
method5 <- qrod_mle(y = data2_outlier[, 1], x = x_model, tau = 0.5, error = 1e-06,
                    iter = 1000, method = "qfunction")
method6 <- qrod_bayes(y = data2_outlier[, 1], x = as.matrix(data2_outlier[, -1]),
                      tau = 0.5, M = 1000, method = "bayes.prob")

p4 <- ggplot(method4, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.5")

p5 <- ggplot(method5, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.5")

p6 <- ggplot(method6, aes(x = case, y = result)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.5")

method7 <- qrod_mle(y = data2_outlier[, 1], x = x_model, tau = 0.9, error = 1e-06,
                    iter = 1000, method = "cook.distance")
method8 <- qrod_mle(y = data2_outlier[, 1], x = x_model, tau = 0.9, error = 1e-06,
                    iter = 1000, method = "qfunction")
method9 <- qrod_bayes(y = data2_outlier[, 1], x = as.matrix(data2_outlier[, -1]),
                      tau = 0.9, M = 1000, method = "bayes.prob")
p7 <- ggplot(method7, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.9")

p8 <- ggplot(method8, aes(x = case, y = distance)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.9")

p9 <- ggplot(method9, aes(x = case, y = result)) +
  geom_point(colour = "darkblue") +
  ggtitle("tau = 0.9")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)






































