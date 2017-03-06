## ---- library
library(quokar)
library(quantreg)
library(ggplot2)
library(gridExtra)
library(purrr)
library(tidyr)
library(dplyr)
library(rggobi)

##---- simplex_method
data(ais)
tau <- c(0.1, 0.5, 0.9)
ais_female <- ais[103:202, ]
br <- rq(BMI ~ LBM, tau = tau, data = ais_female, method = 'br')
coef <- br$coef
br_result <- frame_br(br, tau)
origin_obs <- br_result$all_observation
use_obs <- br_result$fitting_point
ggplot(origin_obs,
    aes(x = value, y = y)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = coef[2, 1], intercept = coef[1,1],
                colour = "gray") +
    geom_abline(slope = coef[2, 2], intercept = coef[1,2],
                colour = "gray") +
    geom_abline(slope = coef[2, 3], intercept = coef[1,3],
                colour = "grey") +
    ylab('y') +
    xlab('x') +
    facet_wrap(~variable, scales = "free_x", ncol = 2) +
    geom_point(data = use_obs, aes(x = value, y = y,
                                        group = tau_flag,
                                        colour = tau_flag,
                                        shape = obs))

br <- rq(BMI ~ LBM + Ht, tau = tau, data = ais_female, method = 'br')
coef <- br$coef
br_result <- frame_br(br, tau)
origin_obs <- br_result$all_observation
use_obs <- br_result$fitting_point
ggplot(origin_obs,
       aes(x = value, y = y)) +
  geom_point(alpha = 0.1) +
  ylab('y') +
  xlab('x') +
  facet_wrap(~variable, scales = "free_x", ncol = 2) +
  geom_point(data = use_obs, aes(x = value, y = y,
                                 group = tau_flag,
                                 colour = tau_flag,
                                 shape = obs))
##---- generate_data
simvar <- function(x, n = 10, method = "grid") UseMethod("simvar")
simvar.factor <- function(x, n = 10, method = "grid"){
  switch(method,
         random = x[sample(length(x), n, replace = TRUE)],
         factor(levels(x), levels = levels(x))
  )
}

simvar.numeric <- function(x, n = 10, method = "grid"){
  rng <- range(x)
  switch(method,
         random = runif(n, rng[1], rng[2]),
         seq(rng[1], rng[2], length = n))
}

generate_data <- function(data, n = 1000, method = "grid"){
  if(method != "random"){
    n <- floor(n ^ (1/ncol(data)))
    df <- data.frame(expand.grid(lapply(data, simvar, n = n,
                                        method = "grid")))
    if(method == "nonaligned"){
      cont <- !sapply(df, is.factor)
      ranges <- lapply(df[, cont], function(x) diff(range(x)))
      df[,cont] <- df[,cont] +
        do.call(cbind, lapply(ranges, function(rng)
          runif(-rng/(2*n), rng/(2*n), n=nrow(df))))
    }
    df
  }else{
    data.frame(sapply(data, simvar, n=n, method=method))
  }
}

##---- multi-variable_br
br_ggobi <- function(data, model, tau){
  n <- nrow(data)
  idx_y = which(colnames(data) == all.vars(model)[1])
  idx_x = which(colnames(data) %in% all.vars(model)[-1])
  y <- as.matrix(data[idx_y])
  x <- as.matrix(data[idx_x])
  object <- rq(model, tau, method = 'br', data = data)
  ntau <- length(tau)
  br_flag <- rep("non-use", n)
  for(i in 1:ntau){
    points <- frame_br(object, tau[i])$fitting_point$index
    br_flag[points] <- paste("use", tau[i], sep="")
  }
  br_data <- cbind(data, br_flag)
}
actual_data <- ais_female[, c('BMI', 'LBM', 'Ht')]
actual_data$flag <- "observations"
taus <- c(0.1, 0.5, 0.9)
model <- BMI ~ LBM + Ht
data <- ais_female
for(tau in taus){
  br_flag <- paste0("use",tau)
  br <- rq(BMI ~ LBM + Ht, tau = tau, data = ais_female,
           method = "br")
  br_data <- br_ggobi(data, model, tau)
  equation_data <- br_data[which(br_data$br_flag == br_flag), ]
  left_side <- data.frame(LBM = equation_data$LBM,
                          Ht = equation_data$Ht)
  br_low <- generate_data(left_side, n = 1000, method = "random")
  right_side <- as.matrix(cbind(1, br_low)) %*% br$coef
  sim_data <- cbind(right_side,br_low)
  sim_data$flag <- paste0("plane",tau)
  colnames(sim_data) <- colnames(actual_data)
  actual_data <- plyr::rbind.fill(actual_data, sim_data)
  actual_data$flag[which(br_data$br_flag == br_flag)] <-
      paste0("plane",tau)
}
g <- ggobi(actual_data)
d <- g[1]
glyph_color(d) <- as.numeric(as.factor(actual_data$flag))

##---- fn_method
tau <- c(0.1, 0.5, 0.9)
fn <- rq(BMI ~ LBM + Ht, data = ais_female, tau = tau, method = 'fn')
fn_obs <- frame_fn_obs(fn, tau)
##For tau = 0.1, plot the observations used in quantile regression
##fitting based on interior point method
fn1 <- fn_obs[,1]
case <- 1: length(fn1)
fn1 <- cbind(case, fn1)
m <- data.frame(y = ais_female$BMI, x1 = ais_female$LBM,
                x2 =ais_female$Ht, fn1)
p <- length(attr(fn$coefficients, "dimnames")[[1]])
m_f <- m %>% gather(variable, value, -case, -fn1, -y)
mf_a <- m_f %>%
  group_by(variable) %>%
  arrange(variable, desc(fn1)) %>%
  filter(row_number() %in% 1:p)
p1 <- ggplot(m_f, aes(x = value, y = y)) +
 geom_point(alpha = 0.1) +
  geom_point(data = mf_a, size = 3) +
  facet_wrap(~variable, scale = "free_x") +
  xlab("x")
 ## For tau = 0.5, plot the observations used in quantile regression
 ##fitting based on interior point method
 fn2 <- fn_obs[,2]
 case <- 1: length(fn2)
 fn2 <- cbind(case, fn2)
 m <- data.frame(y = ais_female$BMI, x1 = ais_female$LBM,
                 x2 = ais_female$Ht, fn2)
 p <- length(attr(fn$coefficients, "dimnames")[[1]])
 m_f <- m %>% gather(variable, value, -case, -fn2, -y)
 mf_a <- m_f %>%
    group_by(variable) %>%
    arrange(variable, desc(fn2)) %>%
    filter(row_number() %in% 1:p )
 p2 <- ggplot(m_f, aes(x = value, y = y)) +
    geom_point(alpha = 0.1) +
    geom_point(data = mf_a, size = 3, colour = "blue", alpha = 0.5) +
    facet_wrap(~variable, scale = "free_x") +
    xlab("x")
 ## For tau = 0.9
 fn3 <- fn_obs[ ,3]
 case <- 1: length(fn3)
 fn3 <- cbind(case, fn3)
 m <- data.frame(y = ais_female$BMI, x1 = ais_female$LBM,
                 x2 = ais_female$Ht, fn3)
 p <- length(attr(fn$coefficients, "dimnames")[[1]])
 m_f <- m %>% gather(variable, value, -case, -fn3, -y)
 mf_a <- m_f %>%
   group_by(variable) %>%
   arrange(variable, desc(fn3)) %>%
   filter(row_number() %in% 1:p )
 p3 <- ggplot(m_f, aes(x = value, y = y)) +
   geom_point(alpha = 0.1) +
   geom_point(data = mf_a, size = 3) +
   facet_wrap(~variable, scale = "free_x") +
   xlab("x")
 grid.arrange(p1, p2, p3, ncol = 1)

##display the weighting matrix
 obs <- data.frame(cbind(fn_obs,id = 1:nrow(fn_obs)))
 selected <- NULL
 for(i in 1:3){
   data <- obs[order(obs[,i],decreasing = T),c(i,4)][1:3,]
   data <- cbind(data,idx=1:3)
   colnames(data) <- c("value","id","idx")
   data = cbind(data,type=rep(colnames(obs)[i],3))
   if(is.null(selected)){
     selected = data
   }else{
     selected =  rbind(selected,data)
   }
 }
 selected$value = round(selected$value,3)
 ggplot(selected,aes(x=idx,y=value,colour=type))+
   geom_point(aes(size=value),alpha=0.5)+
   geom_text(aes(label = id), position=
               position_jitter(width = 0.04, height = 0.04))+
   facet_wrap( ~ type,scale="free_y")
##############################################################
##---- move_y
x <- sort(runif(100))
y <- 40*x + x*rnorm(100, 0, 10)
#locate the outliers
selectedX <- sample(50:100,5)
y2<- y
y2[selectedX] <- x[1:5]*rnorm(5, 0, 10)
y3 <- y2
y3[selectedX] <- y2[selectedX] - 10
y4 <- y3
y4[selectedX] <- y3[selectedX] - 10
df <- data.frame(x, y, y2, y3, y4)
df_m <- df %>% gather(variable, value, -x)
ggplot(df_m, aes(x = x, y=value)) +
  geom_point() +
  xlab("x") +
  ylab("y") +
  facet_wrap(~variable, ncol=2, scale = "free_y") +
  geom_quantile(quantiles = seq(0.1, 0.9, 0.1))

coefs <- 2:5 %>%
  map(~ rq(df[, .] ~ x, data = df, seq(0.1, 0.9, 0.1))) %>%
  map_df(~ as.data.frame(t(as.matrix(coef(.)))))
colnames(coefs) <- c("intercept", "slope")
tau <- rep(seq(0.1, 0.9, by = 0.1), 4)
model <- paste('rq', rep(1:4, each = 9), sep="")

df_m1 <- data.frame(model, tau, coefs)
df_mf <- df_m1 %>% gather(variable, value, -c(model, tau))
ggplot(df_mf, aes(x = tau, y = value, colour = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y") +
  xlab('quantiles') +
  ylab('coefficients')

##---- lm_case
coef_lm <- 2:5 %>%
  map(~ lm(df[, .] ~ x, data = df)) %>%
  map_df(~ as.data.frame(t(as.matrix(coef(.)))))
colnames(coef_lm) <- c("intercept", "slope")

ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = coef_lm$intercept[1], slope = coef_lm$slope[1],
              colour = "red") +
  geom_abline(intercept = coef_lm$intercept[2], slope = coef_lm$slope[2],
              colour = "blue") +
  geom_abline(intercept = coef_lm$intercept[3], slope = coef_lm$slope[3],
              colour = "orange") +
  geom_abline(intercept = coef_lm$intercept[4], slope = coef_lm$slope[4],
              colour = "green")

##---- move_x
#locate the outliers
selectedIdx <- sample(50:100,5)
df <- data.frame(y)
df$y2 <- y
df$x <- x
df$y2[selectedIdx] <- df$x[1:5]*rnorm(5, 0, 10)
df$x2 <- x
df$x2[selectedIdx] <- df$x2[selectedIdx] + 0.2
df$x3 <- df$x2
df$x3[selectedIdx] <- df$x3[selectedIdx] + 0.2
df$x4 <- df$x3
df$x4[selectedIdx] <- df$x4[selectedIdx] + 0.2

df_m <- df %>% gather(variable, value, -y, -y2)
ggplot(df_m, aes(x = value, y=y2)) +
    geom_point() +
    xlab("x") +
    ylab("y") +
    facet_wrap(~variable, ncol=2, scale = "free_x") +
    geom_quantile(quantiles = seq(0.1, 0.9, 0.1))

coefs <- 3:6 %>%
       map(~ rq(df$y2 ~ df[, .], data = df, seq(0.1, 0.9, 0.1))) %>%
       map_df(~ as.data.frame(t(as.matrix(coef(.)))))
colnames(coefs) <- c("intercept", "slope")

tau <- rep(seq(0.1, 0.9, by = 0.1), 4)
model <- paste('rq', rep(1:4, each = 9), sep="")
df_m1 <- data.frame(model, tau, coefs)
df_mf <- df_m1 %>% gather(variable, value, -c(model, tau))
ggplot(df_mf, aes(x = tau, y = value, colour = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y") +
  xlab('quantiles') +
  ylab('coefficients')
###################################################
##---- o-in-number
selectedX1 <- sample(50:100, 5)
y_number_y1 <- y
y_number_y1[selectedX1] <- x[1:5]*rnorm(5, 0, 10)
selectedX2 <- sample(50:100, 10)
y_number_y2 <- y
y_number_y2[selectedX2] <- x[1:10]*rnorm(10, 0, 10)
selectedX3 <- sample(50:100, 15)
y_number_y3 <- y
y_number_y3[selectedX3] <- x[1:15]*rnorm(15, 0, 10)
df <- data.frame(x, y, y_number_y1, y_number_y2, y_number_y3)
df_m <- df %>% gather(variable, value, -x)
ggplot(df_m, aes(x=x, y=value)) +
    geom_point() +
    xlab("x") +
    ylab("y") +
    facet_wrap(~variable, ncol=2) +
    geom_quantile(quantiles = seq(0.1, 0.9, 0.1))
coefs <- 2:5 %>%
            map(~ rq(df[, .] ~ x, data = df, seq(0.1, 0.9, 0.1))) %>%
            map_df(~ as.data.frame(t(as.matrix(coef(.)))))
tau <- rep(seq(0.1, 0.9, by = 0.1), 4)
model <- paste('rq', rep(1:4, each = 9), sep="")
df_m1 <- data.frame(cbind(model, tau, coefs))
ggplot(df_m1, aes(x = tau, y = x, colour = model)) +
    geom_point() +
    geom_line() +
    xlab('quantile') +
    ylab('coefficients')

##---- ALD
data(ais)
x <- matrix(ais$LBM, ncol = 1)
y <- ais$BMI
tau = c(0.1, 0.5, 0.9)
ald_data <- frame_ald(y, x, tau, smooth = 10, error = 1e-6,
                   iter = 2000)
ggplot(ald_data) +
    geom_line(aes(x = r, y = d, group = obs, colour = tau_flag)) +
    facet_wrap(~tau_flag, ncol = 1) +
    xlab('') +
    ylab('Asymmetric Laplace Distribution Density Function')

##---- RR


##---- GCD


##---- QD


##---- BP


##---- BKL

