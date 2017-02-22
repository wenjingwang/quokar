###====== Nonlinear quantile regression with an interior point algorithm ======
## This code is in package `quantreg`
## We modified it to add the weighting matrix of non-linear quantile regression
## process.
"nlrq.control" <- function (maxiter=100, k=2, InitialStepSize = 1,
	big=1e+20, eps=1.0e-07, beta=0.97)
{
    list(maxiter=maxiter, k=k, InitialStepSize = InitialStepSize,
	big=big, eps=eps, beta=beta)
}

"nlrqModel" <- function (form, tau, data, start)
{
    thisEnv <- environment()
    env <- new.env(parent = environment(form))
    for (i in names(data)) {
        assign(i, data[[i]], envir = env)
    }
    ind <- as.list(start)
    parLength <- 0
    for (i in names(ind)) {
        temp <- start[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp, envir = env)
        ind[[i]] <- parLength + seq(along = start[[i]])
        parLength <- parLength + length(start[[i]])
    }
    useParams <- rep(TRUE, parLength)
    lhs <- eval(form[[2]], envir = env)
    rhs <- eval(form[[3]], envir = env)
    resid <- lhs - rhs
    tau <- tau
    dev <- sum(tau * pmax(resid, 0) + (tau - 1) * pmin(resid, 0))
    if (is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function() numericDeriv(form[[3]],
            names(ind), env)
        getRHS <- getRHS.noVarying
        rhs <- getRHS()
    }
    else {
        getRHS.noVarying <- function() eval(form[[3]], envir = env)
        getRHS <- getRHS.noVarying
    }
    dimGrad <- dim(attr(rhs, "gradient"))
    marg <- length(dimGrad)
    if (marg > 0) {
        gradSetArgs <- vector("list", marg + 1)
        for (i in 2:marg) gradSetArgs[[i]] <- rep(TRUE, dimGrad[i -
            1])
        useParams <- rep(TRUE, dimGrad[marg])
    }
    else {
        gradSetArgs <- vector("list", 2)
        useParams <- rep(TRUE, length(attr(rhs, "gradient")))
    }
    npar <- length(useParams)
    gradSetArgs[[1]] <- (~attr(ans, "gradient"))[[2]]
    gradCall <- switch(length(gradSetArgs) - 1, call("[", gradSetArgs[[1]],
        gradSetArgs[[2]]), call("[", gradSetArgs[[1]], gradSetArgs[[2]],
        gradSetArgs[[2]]), call("[", gradSetArgs[[1]], gradSetArgs[[2]],
        gradSetArgs[[2]], gradSetArgs[[3]]), call("[", gradSetArgs[[1]],
        gradSetArgs[[2]], gradSetArgs[[2]], gradSetArgs[[3]],
        gradSetArgs[[4]]))
    getRHS.varying <- function() {
        ans <- getRHS.noVarying()
        attr(ans, "gradient") <- eval(gradCall)
        ans
    }
    QR <- qr(attr(rhs, "gradient"))
    qrDim <- min(dim(QR$qr))
    if (QR$rank < qrDim)
        stop("singular gradient matrix at initial parameter estimates")
    getPars.noVarying <- function() unlist(setNames(lapply(names(ind),
        get, envir = env), names(ind)))
    getPars.varying <- function() unlist(setNames(lapply(names(ind),
        get, envir = env), names(ind)))[useParams]
    getPars <- getPars.noVarying
    internalPars <- getPars()
    setPars.noVarying <- function(newPars) {
        assign("internalPars", newPars, envir = thisEnv)
        for (i in names(ind)) {
            assign(i, unname(newPars[ind[[i]]]), envir = env)
        }
    }
    setPars.varying <- function(newPars) {
        internalPars[useParams] <- newPars
        for (i in names(ind)) {
            assign(i, unname(internalPars[ind[[i]]]), envir = env)
        }
    }
    setPars <- setPars.noVarying
    on.exit(remove(i, data, parLength, start, temp))
    m <- list(resid = function() resid, fitted = function() rhs,
        formula = function() form, tau = function() tau, deviance = function() dev,
        gradient = function() attr(rhs, "gradient"), incr = function() qr.coef(QR, resid), setVarying = function(vary = rep(TRUE,
            length(useParams))) {
            assign("useParams", if (is.character(vary)) {
                temp <- logical(length(useParams))
                temp[unlist(ind[vary])] <- TRUE
                temp
            } else if (is.logical(vary) && length(vary) != length(useParams)) stop("setVarying : vary length must match length of parameters") else {
                vary
            }, envir = thisEnv)
            gradCall[[length(gradCall)]] <<- useParams
            if (all(useParams)) {
                assign("setPars", setPars.noVarying, envir = thisEnv)
                assign("getPars", getPars.noVarying, envir = thisEnv)
                assign("getRHS", getRHS.noVarying, envir = thisEnv)
                assign("npar", length(useParams), envir = thisEnv)
            } else {
                assign("setPars", setPars.varying, envir = thisEnv)
                assign("getPars", getPars.varying, envir = thisEnv)
                assign("getRHS", getRHS.varying, envir = thisEnv)
                assign("npar", length((1:length(useParams))[useParams]),
                  envir = thisEnv)
            }
        }, changeTau = function(newTau) {
            assign("tau", newTau, envir = thisEnv)
            assign("dev", sum(tau * pmax(resid, 0) + (tau - 1) * pmin(resid, 0)), envir = thisEnv)
            return(dev)
        }, setPars = function(newPars) {
            setPars(newPars)
            assign("resid", lhs - assign("rhs", getRHS(), envir = thisEnv),
                envir = thisEnv)
            assign("dev", sum(tau * pmax(resid, 0) + (tau - 1) * pmin(resid, 0)), envir = thisEnv)
            assign("QR", qr(attr(rhs, "gradient")), envir = thisEnv)
            return(QR$rank < min(dim(QR$qr)))
        }, getPars = function() getPars(), getAllPars = function() getPars(),
        getEnv = function() env, trace = function() cat(format(dev),
            ": ", format(getPars()), "\n"), Rmat = function() qr.R(QR),
        predict = function(newdata = list(), qr = FALSE) {
            Env <- new.env()
            for (i in objects(envir = env)) {
                assign(i, get(i, envir = env), envir = Env)
            }
            newdata <- as.list(newdata)
            for (i in names(newdata)) {
                assign(i, newdata[[i]], envir = Env)
            }
            eval(form[[3]], envir = Env)
        })
    class(m) <- "nlrqModel"
    m
}

nlrq_m <- function (formula, data=parent.frame(), start, tau=0.5,
	control, trace=FALSE, method = "L-BFGS-B")
{
    mf <- match.call()
    formula <- as.formula(formula)
    varNames <- all.vars(formula)
    if (length(formula) == 2) {
        formula[[3]] <- formula[[2]]
        formula[[2]] <- 0
    }
    if (missing(start)) {
        if (!is.null(attr(data, "parameters"))) {
            pnames <- names(attr(data, "parameters"))
        }
        else {
            cll <- formula[[length(formula)]]
            func <- get(as.character(cll[[1]]))
            pnames <- as.character(as.list(match.call(func, call = cll))[-1][attr(func, "pnames")])
        }
    }
    else {
        pnames <- names(start)
    }
    varNames <- varNames[is.na(match(varNames, pnames, nomatch = NA))]
    varIndex <- sapply(varNames, function(varName, data, respLength) {
        length(eval(as.name(varName), data))%%respLength == 0
    }, data, length(eval(formula[[2]], data)))
    mf$formula <- parse(text = paste("~", paste(varNames[varIndex], collapse = "+")))[[1]]
    mf$start <- mf$tau <- mf$control <- mf$algorithm <- mf$trace <- mf$method <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- as.list(eval(mf, parent.frame()))
    if (missing(start)) {
        start <- getInitial(formula, mf)
    }
    for (var in varNames[!varIndex]) mf[[var]] <- eval(as.name(var), data)
    ctrl <- nlrq.control()
    if (!missing(control)) {
        control <- as.list(control)
        ctrl[names(control)] <- control
    }
    m <- nlrqModel(formula, tau, mf, start)
    nlrq.calc <- function (model, ctrl, trace) {
        meketon <- function(x, y, w, tau, ctrl) {
            yw <- ctrl$big
            k <- 1
            while(k <= ctrl$k & yw - crossprod(y, w) > ctrl$eps) {
                d <- pmin(tau - w, 1 - tau + w)
                z <- lsfit(x, y, d^2, intercept=FALSE)
                yw <- sum(tau * pmax(z$resid, 0) + (tau - 1) * pmin(z$resid, 0))
                k <- k + 1
                s <- z$resid * d^2
                alpha <- max(ctrl$eps, pmax(s/(tau - w), -s/(1 - tau + w)))
                w <- w + (ctrl$beta/alpha) * s
            }
            coef <- z$coef
            return(list(coef=coef, w=w, d = d))
        }
        model.step <- function(lambda, Step, model, pars) {
            model$setPars(pars + lambda * Step)
            model$deviance()
        }
        w <- rep(0, length(model$resid()))
        snew <- model$deviance()
        sold <- ctrl$big
        nit <- 0
        if (trace) {
            model$trace()
            optim.ctrl <- list(trace=1)
        } else {
            optim.ctrl <- list(trace=0)
        }
        lam0 <- ctrl$InitialStepSize
        D <- list()
        while(sold - snew > ctrl$eps & nit < ctrl$maxiter) {
            z <- meketon(model$gradient(),as.vector(model$resid()), w, tau=tau, ctrl=ctrl)
            Step <- z$coef
            Pars <- model$getPars()
            lam <- try(optim(par=lam0, fn=model.step, method=method, lower=0, upper=1,
	        Step=Step, model=model, pars=Pars, control=optim.ctrl)$par)
            if(inherits(lam,"try.error") || !is.finite(lam))
                stop("optim unable to find valid step size")
            if (trace) {cat("lambda =", lam, "\n")}
            model$setPars(Pars + lam * Step)
            sold <- snew
            snew <- model$deviance()
            w <- qr.resid(qr(model$gradient()), z$w)
            w1 <- max(pmax(w, 0))
            if(w1 > tau) {w <- (w * tau)/(w1 + ctrl$eps)}
            w0 <- max(pmax( - w, 0))
            if(w0 > 1 - tau) {w <- (w * (1 - tau))/(w0 + ctrl$eps)}
            if (trace) {model$trace()}
            if (R.Version()$os == "Win32") {flush.console()}
            nit <- nit + 1
             D[[nit]] <- z$d
        }
        Rho <- function(u,tau) u * (tau - (u < 0))
    	  model$rho <- sum(Rho(model$resid(),tau))
        model$D <- D
        model
    }
    nlrq.out <- list(m=nlrq.calc(m, ctrl, trace), data=substitute(data),
	call=match.call(), PACKAGE = "quantreg")
    nlrq.out$call$control <- ctrl
    nlrq.out$call$trace <- trace
    class(nlrq.out) <- "nlrq"
    nlrq.out
}
