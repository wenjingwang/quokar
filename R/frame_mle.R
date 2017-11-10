globalVariables(c("variable", "value"))
#' @title General Cook's distance or Q-function distance of quantile regression
#' @description dataframe used to plot generalized Cook's distance or
#' Q-function distance for observations.
#' @param y vector, dependent variable for quantile regression
#' @param x matrix, design matrix for quantile regression. For quantile regression
#' model with intercept, the firt column of x is 1.
#' @param tau sigular or vector, quantiles
#' @param error the EM algorithm accuracy of error used in MLE estimation
#' @param iter the iteration frequancy for EM algorithm used in MLE estimation
#' @param method use method 'cook.distance' or 'qfunction'
#' @return generalized Cook's distance or Q-function distance for
#' multiple quantiles
#' @author Wenjing Wang<wenjingwangr@gmail.com>
#' @details
#' Gerneralized Cook's distance and Q-function distance are commonly used in
#' detecting the influence data point when performing regression
#' analysis. They involve the log-likelihood function and estimations of
#' based on the complete and case-deletion data. We used EM algorithm to
#' estimate the coefficiences of quantile regression with asymmetric Laplace
#' distribution.
#' @export
#' @examples
#' library(ggplot2)
#' data(ais)
#' ais_female <- subset(ais, Sex == 1)
#' y <- ais_female$BMI
#' x <- cbind(1, ais_female$LBM, ais_female$Bfat)
#' tau <- c(0.1, 0.5, 0.9)
#' case <- rep(1:length(y), length(tau))
#' GCD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
#'                  method = 'cook.distance')
#' GCD_m <- cbind(case, GCD)
#' ggplot(GCD_m, aes(x = case, y = value )) +
#'   geom_point() +
#'   facet_wrap(~variable, scale = 'free') +
#'   geom_text(data = subset(GCD_m, value > mean(value) + 2*sd(value)),
#'             aes(label = case), hjust = 0, vjust = 0) +
#'   xlab("case number") +
#'   ylab("Generalized Cook Distance")
#'
#' QD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
#'                 method = 'qfunction')
#' QD_m <- cbind(case, QD)
#' ggplot(QD_m, aes(x = case, y = value)) +
#'   geom_point() +
#'   facet_wrap(~variable, scale = 'free')+
#'   geom_text(data = subset(QD_m, value > mean(value) + sd(value)),
#'             aes(label = case), hjust = 0, vjust = 0) +
#'   xlab('case number') +
#'   ylab('Qfunction Distance')
#'

frame_mle <- function(y, x, tau, error = 1e-06,
                      iter = 100, method = c("cook.distance",
                                             "qfunction")){
    if(!is.vector(y)){
      stop("y should be vector")
    }
    if(!is.matrix(x)){
      stop("x should be matrix")
    }
    method <- match.arg(method)
    ntau <- length(tau)
    if(method == 'cook.distance'){
    distances <- qrod_mle(y, x, tau, error,
                          iter, method = "cook.distance")
    }else if(method == 'qfunction'){
      distances <- qrod_mle(y, x, tau, error,
                            iter, method = "qfunction")
    }
    distance_m <- matrix(Reduce(c, distances),
                     ncol = ntau, byrow = FALSE)
    distance_f <- data.frame(distance_m)
    colnames(distance_f) <- paste('tau=', tau, sep = "")
    distance_g <- gather(distance_f, variable, value)
}

