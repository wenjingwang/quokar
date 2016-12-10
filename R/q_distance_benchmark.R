#Q Distance bench-mark
library(dplyr)

q_benchmark <- function(Hessian, q, n){
  if(!is.matrix(Hessian)){
    stop(paste("Hessian must be a matix"))
  }
  if(q < 0 | q > 1){
    stop(paste("q must between 0 to 1"))
  }
  lambda <- eigen(Hessian)$values
  lambda_c <- lambda / sum(lambda)
  e <- eigen(Hessian)$vectors
  mid_data <- data.frame(cbind(lambda_c, t(e)))
  mid_data_d <- dplyr::arrange(mid_data, desc(lambda))
  K = sum(sort(abs(lambda) >= q/sqrt(n)))
  vec_k <- mid_data_d[1:K, -1]^2
  lambda_k <- lambda[1:K]
  M <- apply(lambda_k * vec_k, 2, sum)
  mean_M <- sum(lambda_K)
  sd_M <- sd(M)
  critical_value <- c(2 * mean_M, mean_M + 2 * sd_M)
  return(critical_value)
}

Hessian = matrix(rep(25, 25), nrow = 5)
q_benchmark(Hessian, q = 2, 2)

