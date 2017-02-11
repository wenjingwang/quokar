wh <- function(object, tau){
  x <- object$x
  y <- object$y
  tol <- .Machine$double.eps^(2/3)
  eps <- tol
  big <- .Machine$double.xmax
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  ny <- NCOL(y)
  nsol <- 2
  ndsol <- 2
  if (qr(x)$rank < p)
    stop("Singular design matrix")
  if (tau < 0 || tau > 1) {
    nsol <- 3 * n
    ndsol <- 3 * n
    lci1 <- FALSE
    qn <- rep(0, p)
    cutoff <- 0
    tau <- -1
  }
  else {
    if (p == 1)
      ci <- FALSE
    if (ci) {
      lci1 <- TRUE
      if (tcrit)
        cutoff <- qt(1 - alpha/2, n - p)
      else cutoff <- qnorm(1 - alpha/2)
      if (!iid) {
        h <- bandwidth.rq(tau, n, hs = TRUE)
        bhi <- rq.fit.br(x, y, tau + h, ci = FALSE)
        bhi <- coefficients(bhi)
        blo <- rq.fit.br(x, y, tau - h, ci = FALSE)
        blo <- coefficients(blo)
        dyhat <- x %*% (bhi - blo)
        if (any(dyhat <= 0)) {
          pfis <- (100 * sum(dyhat <= 0))/n
          warning(paste(pfis, "percent fis <=0"))
        }
        f <- pmax(eps, (2 * h)/(dyhat - eps))
        qn <- rep(0, p)
        for (j in 1:p) {
          qnj <- lm(x[, j] ~ x[, -j] - 1, weights = f)$resid
          qn[j] <- sum(qnj * qnj)
        }
      }
      else qn <- 1/diag(solve(crossprod(x)))
    }
    else {
      lci1 <- FALSE
      qn <- rep(0, p)
      cutoff <- 0
    }
  }
  z <- .Fortran("rqbr",
                n = as.integer(n),
                p = as.integer(p),
                n5 = as.integer( n + 5),
                p3 = as.integer(p + 3),
                p4 = as.integer(p + 4),
                a = as.double(x),
                b = as.double(y),
                t = as.double(tau),
                tole = as.double(tol),
                flag = as.integer(1),
                coef = double(p),
                resid = double(n),
                s =  integer(n),
                wa = double((n + 5) * (p + 4)),
                wb = double(n),
                nsol = as.integer(nsol),
                ndsol = as.integer(ndsol),
                sol = double((p + 3) * nsol),
                dsol = double(n * ndsol),
                lsol = as.integer(0),
                wh = integer(p * nsol),
                qn = as.double(qn),
                cutoff = as.double(cutoff),
                ci = double(4 * p),
                tnmat = double(4 *  p),
                big =  as.double(big),
                lci1 = as.logical(lci1),
                PACKAGE = "quokar")
  return(z$wh[1:NCOL(x)])
}
