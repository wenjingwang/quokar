#' @title General cook distance of quantile regression
#' @param y dependent variable for quantile regression
#' @param x design matrix for quantile regression
#' @param tau quantiles
#' @param error the EM algorithm accuracy of error used in MLE
#' estimation
#' @param iter the iteration frequancy for EM algorithm used in
#' estimation
#' @return generalized cook distance for multiple quantiles
#' @author Wenjing Wang<wenjingwang1990@gmail.com>
#' @details This function used to prepare the data frame to
#' plot general cook distance.
#' @seealso see documentation of qrod_mle
#' @examples
#' data(ais)
#' y <- ais$BMI
#' x <- cbind(1, ais$LBM))
#' case <- rep(1:length(y), length(tau))
#' GCD <- frame_GCD(y, x, tau, error = 1e-06, iter = 100)
#' GCD_m <- cbind(case, GCD)
#' ggplot(GCD_m, aes(x = case, y = value )) +
#'   geom_point() +
#'   facet_wrap(~variable, scale = 'free') +
#'   xlab("case number") +
#'   ylab("Generalized Cook Distance")
#'
#'

frame_GCD <- function(y, x, tau, error = 1e-06,
                      iter = 100){
    ntau <- length(tau)
    distances <- qrod_mle(y, x, tau, error,
                          iter, method = "cook.distance")
    distance_m <- matrix(Reduce(c, distances),
                     ncol = ntau, byrow = FALSE)
    distance_f <- data.frame(distance_m)
    colnames(distance_f) <- paste('tau=', tau, sep = "")
    distance_g <- distance_f %>% gather(variable, value)
}
head(distance_g)




