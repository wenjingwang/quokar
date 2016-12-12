#'@title Calculating the benchmark of Studentised residuals
#'of quantile regression
#'
#'@description
#' Benchmark of studentized residuals
#'

stud_resid_benchmark <- function(alpha, p, n){
  stud_resid_benchmark <- qt(1-alpha/2, n - 2*p -1, lower.tail = TRUE)
}
