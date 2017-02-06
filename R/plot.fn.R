#' Create plot of the fitting process of quantile regression
#' using 'fn' algorithm
#'
#' @param y dependent variable
#' @param x independent variable
#' @param tau quantile
#'
#'
#' @author Wenjing Wang \email{wenjingwang1990@gmail.com}
#' @useDynLib quokar
#' @export
#'
loadRQFnbj <- function(){
  if(!is.loaded("rqfnb_path")){
    curWd <- getwd()
    setwd(tempdir())
    cat("## we use the BLAS and now also the LAPACK library:
        PKG_LIBS= $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)",
        file="Makevars")
    system2(command=R.home("bin/R"),
            args=c("CMD","SHLIB","rqfnb_path.f"))
    dyn.load(paste0("rqfnb_path",.Platform$dynlib.ext))
    setwd(curWd)
  }
}
loadRQFnbj()

source("loadRQFnbj.r")

data(ais)
y <- ais$BMI
x <- cbind(1, ais$LBM)
tau <- 0.5
beta <- 0.99995
eps <- 1e-06
n <- length(y)
x <- as.matrix(x)
p <- ncol(x)
if (n != nrow(x))
  stop("x and y don't match n")
if (tau < eps || tau > 1 - eps)
  stop("No parametric Frisch-Newton method.
       Set tau in (0,1)")
rhs <- (1 - tau) * apply(x, 2, sum)
d <- rep(1, n)
u <- rep(1, n)
wn <- rep(0, 10 * n)
wn[1:n] <- (1 - tau)
z <- .Fortran("rqfnbj",
              n = as.integer(n),
              p = as.integer(p),
              a = as.double(t(as.matrix(x))),
              c = as.double(-y),
              rhs = as.double(rhs),
              d = as.double(d),
              u = as.double(u),
              beta = as.double(beta),
              eps = as.double(eps),
              wn = as.double(wn),
              wp = double((p + 3) * p),
              it.count = integer(3),
              info = integer(1),
              it.routine = double(50*p),
              PACKAGE = "quokar"
)
print(-matrix(z$it.routine,ncol=p,byrow=TRUE))
