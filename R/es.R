es <- function(object){

  es_data <- matrix(nrow = length(object$y), ncol = ncol(object$x) - 1)

  for(i in 1: length(object$y)){

    if(object$residual[i] <= sort(object$residual)[ncol(object$x)]){

      es_data[i, ] <- object$x[i, 2:ncol(object$x)]

    }
  }

  return(es_data)
}
