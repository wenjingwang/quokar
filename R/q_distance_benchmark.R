#'@title Calculating the benchmark of q function distance
#'@param Hessian the hessian matrix for benchmark calculating
#'@param q the parameter indicate the location of the local
#'influenctial observations
#'@description  This function calculate the benchmark of local influence assessing based
#'on Q function.
#'@details
#'We propose the following Q-displacement function, because
#'of the intractable likelihood function. The Q-displacement function is:
#'
#'\deqn{f_{Q}(w) = 2[Q(\hat{\theta}|\hat{\theta})-Q{\hat{\theta}(w)|\hat{\theta}}]}
#'
#'where \eqn{\hat{\theta}(w)} is the estimate of \eqn{\theta} which maximizes
#'
#'\deqn{Q(\theta, w|\hat{\theta})=E{L_{c}(\theta, w|Y_{c})|Y_{o}, \hat{\theta}}}
#'
#'The benchmark of the q function distance is \eqn{\={M}(m_{0})+2S_{M}(m_{0})}
#'
#'where,
#'\deqn{-2\"{Q} = \sum{\lambda_{i}e_{i}e_{i}^{'}}}
#'\deqn{M(m_{0}) = \sum {\~{\lambda_{i}}e_{i}^{2}}}
#'
#'@importFrom dplyr arrange desc
#'
#'@references Zhu H T, Lee S Y. Local influence for incomplete data models[J]
#'. Journal of the Royal Statistical Society: Series B (Statistical Methodology)
#', 2001, 63(1): 111-126.
#'
q_benchmark <- function(Hessian, q){
  if(!is.matrix(Hessian)){
    stop("Hessian must be a matix")
  }
  if(q < 0 | q > 1){
    stop("q must larger than 0")
  }
  n <- nrow(Hessian)
  lambda <- eigen(Hessian)$values
  lambda_c <- lambda / sum(lambda)
  e <- eigen(Hessian)$vectors
  mid_data <- data.frame(cbind(lambda_c, t(e)))
  mid_data_d <- arrange(mid_data, desc(lambda))
  k = sum(sort(abs(lambda) >= q/n))
  vec_k <- mid_data_d[1:k, -1]^2
  lambda_k <- lambda[1:k]
  M <- apply(lambda_k * vec_k, 2, sum)
  mean_M <- sum(lambda_k)
  sd_M <- stats::sd(M)
  critical_value <- c(2 * mean_M, mean_M + 2 * sd_M)
  return(critical_value)
}
