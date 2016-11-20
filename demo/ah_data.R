library(ggplot2)
library(dplyr)
library(plyr)
library(RColorBrewer)

#Analyzing the House Price Data
data(ah)
ggplot(ah, aes(SalePrice)) +
  geom_histogram(aes(y = ..density..), colour = "darkblue") +
  geom_density()
#Analyze the Neighborhood
ggplot(ah, aes(factor(Neighborhood))) +
  geom_bar(stat = "count") +
  xlab("Neighborhood") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0))
#Analyze the Sale Price by Neighborhood
Price_Neighbor <- aggregate(ah$SalePrice, list(ah$Neighborhood), mean)
ggplot(Price_Neighbor, aes(x = Group.1, y = x)) +
  geom_bar(stat = "identity") +
  xlab("Neighborhood") +
  ylab("Average Sale Price") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0))
#Sale Price by year
Price_Year <- aggregate(ah$SalePrice, list(ah$Neighborhood,
                                           factor(ah$YrSold_YYYY)), mean)
ggplot(Price_Year, aes(x=Group.1, y=x, colour=Group.2, group=Group.2)) +
  geom_point() +
  geom_line(lwd = 1) +
  xlab("Neighborhood") +
  ylab("Average Sale Price") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0))
#Lot Area
data1 <- subset(ah, LotArea > 50000)
ggplot(ah, aes(x = LotArea, y = SalePrice)) +
  geom_point(colour = "darkblue") +
  geom_text(data = data1,
            aes(x = LotArea, y = SalePrice - 10000),
            size =3, label = data1$Neighborhood)
#Year of Built
ah1 <- subset(ah, YrBuilt != 0)
Price_YrBuilt <- aggregate(ah1$SalePrice, list(ah1$YrBuilt), mean)
year_built <- 2013 - Price_YrBuilt$Group.1 + 1
Price_YrBuilt <- cbind(Price_YrBuilt, year_built)
ggplot(data = Price_YrBuilt, aes(x = year_built, y = x)) +
       geom_point(colour = "darkblue") +
       geom_smooth() +
       xlab("Built Years") +
       ylab("Average Sale Price")
#House Style
HS <- aggregate(ah$SalePrice, list(ah$HouseStyle), mean)
ggplot(HS, aes(x = Group.1, y = x)) +
  geom_bar(stat = "identity") +
  xlab("House Style") +
  ylab("Sale Price") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0))
#########################################
ah1 <- subset(ah, YrBuilt != 0)
year_built <- 2013 - ah1$YrBuilt + 1
ah2 <- cbind(ah1, year_built)
head(ah2)

ah_lm <- lm(log(SalePrice) ~ factor(HouseStyle) + factor(YrSold_YYYY) +
              year_built, data = ah2)

summary(ah_lm)

n <- nrow(ah2)
y <- log(ah2$SalePrice)
housestyle.f <- factor(ah2$HouseStyle)
dum_hs <- model.matrix( ~ housestyle.f )
yrsold.f <- factor(ah2$YrSold_YYYY)
dum_ys <- model.matrix( ~ yrsold.f)
x <- cbind(rep(1, n), year_built, dum_hs[, 2:ncol(dum_hs)], dum_ys[, 2:ncol(dum_ys)])



tau = 0.1
GCD <- ALDqr_GCD(y, x, tau, error = 1e-06, iter = 2000)
n <- length(y)
case <- 1:n
critical_value <- 2*(tau + 1) / n
ah_GCD <- data.frame(case = case, GCD = GCD)
p1 <- ggplot(ah_GCD, aes(x = case, y = GCD)) +
         geom_point(colour = "darkblue") +
         geom_hline(yintercept = critical_value, colour = "red") +
         geom_text(data = subset(ah_GCD, GCD > critical_value),
                   aes(case, GCD + 0.001, label = case))

tau = 0.2
GCD <- ALDqr_GCD(y, x, tau, error = 1e-06, iter = 2000)
n <- length(y)
case <- 1:n
critical_value <- 2*(tau + 1) / n
ah_GCD <- data.frame(case = case, GCD = GCD)
p2 <- ggplot(ah_GCD, aes(x = case, y = GCD)) +
        geom_point(colour = "darkblue") +
        geom_hline(yintercept = critical_value, colour = "red") +
        geom_text(data = subset(ah_GCD, GCD > critical_value),
                  aes(case, GCD + 0.001, label = case))


tau = 0.5
GCD <- ALDqr_GCD(y, x, tau, error = 1e-06, iter = 2000)
n <- length(y)
case <- 1:n
critical_value <- 2*(tau + 1) / n
ah_GCD <- data.frame(case = case, GCD = GCD)
p3 <- ggplot(ah_GCD, aes(x = case, y = GCD)) +
        geom_point(colour = "darkblue") +
        geom_hline(yintercept = critical_value, colour = "red") +
        geom_text(data = subset(ah_GCD, GCD > critical_value),
                  aes(case, GCD + 0.001, label = case))

tau = 0.9
GCD <- ALDqr_GCD(y, x, tau, error = 1e-06, iter = 2000)
n <- length(y)
case <- 1:n
critical_value <- 2*(tau + 1) / n
ah_GCD <- data.frame(case = case, GCD = GCD)
p4 <- ggplot(ah_GCD, aes(x = case, y = GCD)) +
        geom_point(colour = "darkblue") +
        geom_hline(yintercept = critical_value, colour = "red") +
        geom_text(data = subset(ah_GCD, GCD > critical_value),
                  aes(case, GCD + 0.001, label = case))


grid.arrange(p1, p2, p3, p4, ncol = 2)


tau <- c(0.1, 0.2, 0.5, 0.9)
n <- length(y)
case <- 1: n

QD <- ALDqr_QD(y, x, tau[1], error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p1 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 20),
            aes(case, QD+1, label = case), size = 2)

QD <- ALDqr_QD(y, x, tau=0.2, error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p2 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 10),
            aes(case, QD+1, label = case), size = 2)

QD <- ALDqr_QD(y, x, tau[3], error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p3 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 5),
            aes(case, QD+1, label = case), size = 2)

QD <- ALDqr_QD(y, x, tau[4], error = 1e-06, iter = 2000)
data_QD <- data.frame(case = case, QD = QD)
p4 <- ggplot(data_QD, aes(x = case, y = QD)) +
  geom_point(colour = "darkblue") +
  geom_text(data = subset(data_QD, QD > 20),
            aes(case, QD+1, label = case), size = 2)

grid.arrange(p1, p2, p3, p4, ncol = 2)


#######################
#descriptive statistics
data(ah)
colnames(ah)
ah.omitna <- na.omit(ah) #delete the na value
ah.omitna2 <- subset(ah.omitna, YrBuilt != 0) #delete the 0 value standing for na
year_built <- 2013 - ah.omitna2$YrBuilt + 1
ah2 <- cbind(ah.omitna2, year_built)


ah2$YrSold_YYYY <- as.character(ah2$YrSold_YYYY)
ah2$SaleDate <- as.character(ah2$SaleDate)
ah2$MoSold_MM <- as.character(ah2$MoSold_MM)
ah2$YrBuilt <- as.character(ah2$YrBuilt)

#continous variables:
#install.packages("stargazer")
library(stargazer)
ah.num <- ah2[ ,sapply(ah2, is.numeric)]
stargazer(ah.num)







