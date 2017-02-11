globalVariables(c("variable", "value"))
#' @title General cook distance or Qfunction distance of quantile regression
#' @param y dependent variable for quantile regression
#' @param x design matrix for quantile regression
#' @param tau quantiles
#' @param error the EM algorithm accuracy of error used in MLE
#' estimation
#' @param iter the iteration frequancy for EM algorithm used in
#' estimation
#' @param method use method 'cook.distance' or 'qfunction'
#' @description data frame used to plot generalized cook distance or
#' qfunction distance
#' @return generalized cook distance or qfunction distance for
#' multiple quantiles
#' @author Wenjing Wang<wenjingwang1990@gmail.com>
#' @details This function used to prepare the data frame to
#' plot general cook distance or qfunction distance
#' @seealso see documentation of qrod_mle
#' @importFrom tidyr gather
#' @export
#' @examples
#' data(ais)
#' y <- ais$BMI
#' x <- cbind(1, ais$LBM)
#' tau <- c(0.1, 0.5, 0.9)
#' case <- rep(1:length(y), length(tau))
#' GCD <- frame_mle(y, x, tau, error = 1e-06, iter = 100,
#'                  method = 'cook.distance')
#' GCD_m <- cbind(case, GCD)
#' ggplot(GCD_m, aes(x = case, y = value )) +
#'   geom_point() +
#'   facet_wrap(~variable, scale = 'free') +
#'   xlab("case number") +
#'   ylab("Generalized Cook Distance")
#' QD <- frame_mle(y, x, tau, error = 1e-06, iter = 100,
#'                 method = 'qfunction')
#' QD_m <- cbind(case, QD)
#' ggplot(QD_m, aes(x = case, y = value)) +
#'   geom_point() +
#'   facet_wrap(~variable, scale = 'free')+
#'   xlab('case number') +
#'   ylab('Qfunction Distance')
#'

frame_mle <- function(y, x, tau, error = 1e-06,
                      iter = 100, method){
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

