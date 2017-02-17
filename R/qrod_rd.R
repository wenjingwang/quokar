#' @title Robust Distance for Outlier Detection of Quantile Regression
#'     model
#' @param
#' @param
#' @details It is possible that the original data is in $p$
#' dimensional space, but the $h$-point subset that yields the minimum
#' covariance determinant lies in a lower-dimensional hyperplane.
#' Applying the conoical MCD algorithm to such a data set would result
#' in a singular covariance problem(called exact fit in Rousseeuw and
#' Van Driessen(1999)), so that the relevant robust distances cannot be
#' computed. To deal with the singularity problem and provide further
#' leverage point analysis, we implement a generalized MCD algorithm.
#' The algorithm distinguishes in-hyper plane points from off-hyper
#' plane points, and performs MCD leverage point analysis in the
#' dimension-reduced space by projecting all points onto the
#' hyperplane. The generalized MCD algorithm solves that problem by
#' identifying all male observations as off-plane leverage points, and
#' then carries out the leverage point analysis with all the other
#' covariates being centered separately.
#' We use MCD to detect outliers in quantile regression model.
#'
#'
#' @description
#' @author Wenjing Wang
#'
#'
qrod_rd <- function(y, x, tau){
    ##The robust multivariable location and scale estimates computed
    ##with the minimum covariance determinant(MCD) method of Rousseeuw
    ##and Van Driessen(1999)
    T(A) <-
    C(A) <-
        ##Robust distance is
        RD <-
    ##detecting leverage points
    if(RD <= C){
        leverage_point[i] <- 0
    }else{
        leverage_point[i] <- 1
    }
    ##detecting outliers
    if(r <= k*sigma){
        outlier_point[i] <- 0
    }else{
        outlier_point[i] <- 1
    }

}


























