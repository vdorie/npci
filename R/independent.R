setClass("gpci.independent", contains = "gpci",
         slots = list(deviance.0 = "function",
                      K.0 = "matrix", 
                      L.0 = "matrix",
                      deviance.1 = "function",
                      K.1 = "matrix",
                      L.1 = "matrix"))

setMethod("initialize", "gpci.independent",
  function(.Object, y, x, z)
{
  df <- as.data.frame(x)
  df$.y <- y
  trans <- getTransformations(df)
  df <- transform(df, trans$standardize)
  
  x <- as.matrix(df[,names(df) %not_in% ".y"])
  xt.0 <- t(if (ncol(x) > 1) x[z == 0,] else x[z == 0])
  xt.1 <- t(if (ncol(x) > 1) x[z == 1,] else x[z == 1])
  n.0 <- sum(z == 0)
  n.1 <- sum(z == 1)
  
  .Object@data <- list(y = df$.y,
                       z = z,
                       x = x,
                       y.0 = df$.y[z == 0],
                       y.1 = df$.y[z == 1],
                       xt.0 = xt.0,
                       xt.1 = xt.1)
  .Object@K.0 <- matrix(NA_real_, n.0, n.0)
  .Object@L.0 <- matrix(NA_real_, n.0, n.0)
  .Object@K.1 <- matrix(NA_real_, n.1, n.1)
  .Object@L.1 <- matrix(NA_real_, n.1, n.1)
  
  .Object@trans <- trans
  .Object@env <- new.env(parent = baseenv())
  
  devianceEnv <- args2env(baseenv(), object = .Object)
  devianceEnv$transformPars.0 <- function(x) x
  devianceEnv$transformPars.1 <- function(x) x
  
  .Object@deviance.0 <- deviance.0.gpci.independent
  .Object@deviance.1 <- deviance.1.gpci.independent
  environment(.Object@deviance.0) <- devianceEnv
  environment(.Object@deviance.1) <- devianceEnv
  
  validObject(.Object)
  .Object
})

deviance.0.gpci.independent <- function(pars)
{
  n.covPars <- ncol(object@data$x)
  
  pars <- transformPars.0(pars)
  
  sig_f.0_sq <- exp(pars[1L])
  scales.0 <- exp(-0.5 * pars[-1L])
  
  n.0 <- sum(object@data$z == 0)
  .Call(npci:::C_npci_grouped_updateCovMatrix, object@K.0, object@data$xt.0, object@data$xt.0, scales.0, sig_f.0_sq)
  
  cholResult <- .Call(npci:::C_npci_grouped_updateLeftFactor, object@L.0, object@K.0)
  if (cholResult != 0L) return(.Machine$double.xmax * .Machine$double.eps^2)
  
  xr <- solve(object@L.0, rep(1, n.0))
  xtx <- crossprod(xr)
  mu.0 <- solve(xtx, crossprod(xr, solve(object@L.0, object@data$y.0)))[1]
  
  sig_y.0_sq.hat <- crossprod(solve(object@L.0, object@data$y.0 - mu.0))[1] / n.0
  
  result <- n.0 * (log(sig_y.0_sq.hat) + log(2 * pi)) + 2 * sum(log(diag(object@L.0))) + n.0
  attr(result, "beta") <- mu.0
  attr(result, "sig_y_sq") <- sig_y.0_sq.hat
  result
}

deviance.1.gpci.independent <- function(pars)
{
  n.covPars <- ncol(object@data$x)
  
  pars <- transformPars.1(pars)
  
  sig_f.1_sq <- exp(pars[1L])
  scales.1 <- exp(-0.5 * pars[-1L])
  
  n.1 <- sum(object@data$z == 1)
  .Call(npci:::C_npci_grouped_updateCovMatrix, object@K.1, object@data$xt.1, object@data$xt.1, scales.1, sig_f.1_sq)
  
  cholResult <- .Call(npci:::C_npci_grouped_updateLeftFactor, object@L.1, object@K.1)
  if (cholResult != 0L) return(.Machine$double.xmax * .Machine$double.eps^2)

  xr <- solve(object@L.1, rep(1, n.1))
  xtx <- crossprod(xr)
  mu.1 <- solve(xtx, crossprod(xr, solve(object@L.1, object@data$y.1)))[1]

  sig_y.1_sq.hat <- crossprod(solve(object@L.1, object@data$y.1 - mu.1))[1] / n.1

  result <- n.1 * (log(sig_y.1_sq.hat) + log(2 * pi)) + 2 * sum(log(diag(object@L.1))) + n.1
  attr(result, "beta") <- mu.1
  attr(result, "sig_y_sq") <- sig_y.1_sq.hat
  result
}

optimize.gpci.independent <- function(object, n.testPoints = c(50L, 15L), n.modes = 4L, verbose = FALSE, start.only = FALSE, transform.pars = TRUE)
{
  if (is.null(object@env$opt.0)) object@env$opt.0 <- list()
  if (length(object@env$opt.0) >= n.modes && length(object@env$opt.1) >= n.modes) return(invisible(NULL))
  
  n.covPars <- ncol(object@data$x)
  
  devEnv <- environment(object@deviance.0)
  if (is.null(object@env$optStart.0)) {
    n.testPoints.0 <- n.testPoints[1] + n.testPoints[2L] * n.covPars
    testMins  <- c(0, rep(-3, n.covPars))
    testMaxes <- c(10, rep(3, n.covPars))
    
    testPoints <- matrix(runif(n.testPoints.0 * (n.covPars + 1L), testMins, testMaxes), n.testPoints.0, byrow = TRUE)
    
    if (verbose) cat("calculating deviances for group 0\n")
    testDeviances <- t(sapply(seq_len(nrow(testPoints)), function(i) {
      object@deviance.0(testPoints[i,])
    }))
    
    devOrder <- order(testDeviances)
    testPoints    <- testPoints[devOrder,]
    testDeviances <- testDeviances[devOrder]
    
    n.testPoints.0 <- ceiling(0.4 * n.testPoints.0)
    testPoints    <- testPoints[seq_len(n.testPoints.0),]
    testDeviances <- testDeviances[seq_len(n.testPoints.0)]

    
    kfit <- kmeans(testPoints, n.modes, nstart = 5L)
    object@env$optStart.0 <- kfit$centers
  }
  if (is.null(object@env$optStart.1)) {
    n.testPoints.1 <- n.testPoints[1] + n.testPoints[2L] * n.covPars
    testMins  <- c(0, rep(-3, n.covPars))
    testMaxes <- c(10, rep(3, n.covPars))
    
    testPoints <- matrix(runif(n.testPoints.1 * (n.covPars + 1L), testMins, testMaxes), n.testPoints.1, byrow = TRUE)
    
    if (verbose) cat("calculating deviances for group 1\n")
    testDeviances <- t(sapply(seq_len(nrow(testPoints)), function(i) {
      object@deviance.1(testPoints[i,])
    }))
    
    devOrder <- order(testDeviances)
    testPoints    <- testPoints[devOrder,]
    testDeviances <- testDeviances[devOrder]
    
    n.testPoints.1 <- ceiling(0.4 * n.testPoints.1)
    testPoints    <- testPoints[seq_len(n.testPoints.1),]
    testDeviances <- testDeviances[seq_len(n.testPoints.1)]
    
    kfit <- kmeans(testPoints, n.modes, nstart = 5L)
    object@env$optStart.1 <- kfit$centers
  }
  
  if (identical(start.only, TRUE)) return(invisible(NULL))
  
  startTime <- Sys.time()
  if (length(object@env$opt.0) < n.modes) for (i in seq.int(length(object@env$opt.0) + 1L, n.modes)) {
    if (verbose) cat("optimizing iteration ", i, " for group 0\n", sep = "")
    object@env$opt.0[[length(object@env$opt.0) + 1L]] <- optim(object@env$optStart.0[i,], object@deviance.0, method = "BFGS", hessian = TRUE)
  }
  
  if (length(object@env$opt.1) < n.modes) for (i in seq.int(length(object@env$opt.1) + 1L, n.modes)) {
    if (verbose) cat("optimizing iteration ", i, " for group 1\n", sep = "")
    object@env$opt.1[[length(object@env$opt.1) + 1L]] <- optim(object@env$optStart.1[i,], object@deviance.1, method = "BFGS", hessian = TRUE)
  }
  endTime <- Sys.time()
  if (verbose) print(endTime - startTime)
  
  devs <- sapply(object@env$opt.0, function(opt) opt$value)
  object@env$opt.0 <- object@env$opt.0[order(devs)]
  devs <- sapply(object@env$opt.1, function(opt) opt$value)
  object@env$opt.1 <- object@env$opt.1[order(devs)]

  if (!identical(transform.pars, FALSE)) {
    if (is.null(environment(devEnv$transformPars.0)$L)) {
      tryResult <- tryCatch(L <- solve(chol(object@env$opt.0[[1]]$hessian)), error = function(e) e)
      if (is(tryResult, "error")) L <- diag(1, n.covPars + 1L)
      
      x.0 <- object@env$opt.0[[1]]$par
      devEnv$transformPars.0 <- function(x) L %*% x + x.0
      environment(devEnv$transformPars.0) <- args2env(baseenv(), L, x.0)
      for (i in seq_along(object@env$opt.0))
        object@env$opt.0[[i]]$par <- solve(L, object@env$opt.0[[i]]$par - x.0)
    }
    if (is.null(environment(devEnv$transformPars.1)$L)) {
      tryResult <- tryCatch(L <- solve(chol(object@env$opt.1[[1]]$hessian)), error = function(e) e)
      if (is(tryResult, "error")) L <- diag(1, n.covPars + 1L)
      
      x.0 <- object@env$opt.1[[1]]$par
      devEnv$transformPars.1 <- function(x) L %*% x + x.0
      environment(devEnv$transformPars.1) <- args2env(baseenv(), L, x.0)
      for (i in seq_along(object@env$opt.1))
        object@env$opt.1[[i]]$par <- solve(L, object@env$opt.1[[i]]$par - x.0)
    }

  }
  
  invisible(NULL)
} 

sample.gpci.independent <- function(object, n.chains = 4L, n.iter = 300L, n.burn = 150L, verbose = FALSE) {
  startingPoints.0 <- startingPoints.1 <- "random"
  if (!is.null(object@env$opt.0) && length(object@env$opt.0) > 0) {
    startingPoints.0 <- lapply(seq_along(object@env$opt.0), function(i) {
      temp <- object@deviance.0(object@env$opt.0[[i]]$par)
      list(transformedPars = object@env$opt.0[[i]]$par,
           log_sig_y_sq = log(attr(temp, "sig_y_sq")),
           mu = attr(temp, "beta"))
    })
    startingPoints.0 <- rep_len(startingPoints.0, n.chains)
  }
  if (!is.null(object@env$opt.1) && length(object@env$opt.1) > 0) {
    startingPoints.1 <- lapply(seq_along(object@env$opt.1), function(i) {
      temp <- object@deviance.1(object@env$opt.1[[i]]$par)
      list(transformedPars = object@env$opt.1[[i]]$par,
           log_sig_y_sq = log(attr(temp, "sig_y_sq")),
           mu = attr(temp, "beta"))
    })
    startingPoints.1 <- rep_len(startingPoints.1, n.chains)
  }
  
  devEnv <- environment(object@deviance.0)

  fit.0 <- if (!is.null(object@env$sampler.0)) object@env$sampler.0 else NA
  
  x.0 <- if (ncol(object@data$x) > 1) object@data$x[object@data$z == 0,] else as.matrix(object@data$x[object@data$z == 0])
    
  .tempEnv <- NULL
  if (!exists("cpp_object_initializer", envir = .GlobalEnv)) {
    .tempEnv <- npci:::args2env(.GlobalEnv, cpp_object_initializer = Rcpp::cpp_object_initializer)
    attach(.tempEnv)
  }
  object@env$sampler.0 <-
    stan(model_code = independentStanCode,
         data = list(numObservations = length(object@data$y.0),
                     numPredictors   = ncol(object@data$x),
                     x = x.0,
                     y = object@data$y.0,
                     parL = environment(devEnv$transformPars.0)$L,
                     parM = environment(devEnv$transformPars.0)$x.0),
         fit = fit.0,
         init = startingPoints.0,
         chains = n.chains,
         iter = n.iter,
         warmup = n.burn,
         verbose = verbose)
  
  
  fit.1 <- if (!is.null(object@env$sampler.1)) object@env$sampler.1 else object@env$sampler.0
  x.1 <- if (ncol(object@data$x) > 1) object@data$x[object@data$z == 1,] else as.matrix(object@data$x[object@data$z == 1])
  
  object@env$sampler.1 <-
    stan(model_code = independentStanCode,
         data = list(numObservations = length(object@data$y.1),
                     numPredictors   = ncol(object@data$x),
                     x = x.1,
                     y = object@data$y.1,
                     parL = environment(devEnv$transformPars.1)$L,
                     parM = environment(devEnv$transformPars.1)$x.0),
         fit = fit.1,
         init = startingPoints.1,
         chains = n.chains,
         iter = n.iter,
         warmup = n.burn,
         verbose = verbose)
  
  if (!is.null(.tempEnv)) detach(.tempEnv)
  
  invisible(NULL)
}

predict.gpci.independent <- function(object, x, z, error = NULL, pars = "sampler")
{
  x <- as.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- if (!is.null(colnames(object@data$x)))
      colnames(object@data$x)
    else
      names(object@trans$standardize)[names(object@trans$standardize) %not_in% ".y"]
  }
  
  x <- as.matrix(transform(as.data.frame(x), object@trans$standardize))
  z <- rep_len(z, nrow(x))
  
  ind.0 <- which(z == 0)
  ind.1 <- which(z != 0)
  
  xt.0 <- if (length(ind.0) > 0)
    t(if (ncol(x) > 1) x[ind.0,] else as.matrix(x[ind.0]))
  else
    matrix(NA_real_, 0, 0)
  xt.1 <- if (length(ind.1) > 0) 
    t(if (ncol(x) > 1) x[ind.1,] else as.matrix(x[ind.1]))
  else
    matrix(NA_real_, 0, 0)
  
  devEnv <- environment(object@deviance.0)
  if ((is.null(object@env$opt.0) && is.null(object@env$sampler.0)) ||
      (is.null(object@env$opt.1) && is.null(object@env$sampler.1))) optimize.gpci.independent(object)
  
  if (identical(pars, "sampler")) {
    if (length(ind.0) > 0) {
      KJ <- matrix(NA_real_, ncol(object@data$xt.0), ncol(xt.0))
      J  <- matrix(NA_real_, ncol(xt.0), ncol(xt.0))
      if (is.null(object@env$sampler.0)) sample.gpci.independent(object)
      pars <- extract(object@env$sampler.0, object@env$sampler.0@model_pars[seq_len(length(object@env$sampler.0@model_pars) - 1L)])
      
      origPars <- t(sapply(seq_len(nrow(pars$transformedPars)), function(i) devEnv$transformPars.0(pars$transformedPars[i,])))
      
      sig_f.0_sq <- exp(origPars[,1])
      origPars <- as.matrix(origPars[,-1])
      sig_y.0_sq <- exp(pars$log_sig_y_sq)
      mu <- pars$mu
      
      pst.0 <- lapply(seq_len(nrow(origPars)), function(i) {
        scales <- exp(-0.5 * origPars[i,])
        .Call(C_npci_grouped_updateCovMatrix, object@K.0, object@data$xt.0, object@data$xt.0, scales, sig_f.0_sq)
        
        .Call(C_npci_grouped_updateCovMatrix, KJ, object@data$xt.0, xt.0, scales, sig_f.0_sq)
        
        .Call(C_npci_grouped_updateLeftFactor, object@L.0, object@K.0)
        LKJ <- solve(object@L.0, KJ)
        mu.pst <- transform(.y = as.vector(mu[i] + crossprod(LKJ, solve(object@L.0, object@data$y.0 - mu[i]))),
                            trans = object@trans$inverse, simplify = TRUE)
        
        if (is.null(error)) return(mu.pst)
        
        .Call(C_npci_grouped_updateCovMatrix, J, xt.0, xt.0, scales, sig_f.0_sq)
        vcov.pst <- sig_y.0_sq[i] * (J - crossprod(LKJ))
        if (identical(error, "ppd")) vcov.pst <- vcov.pst + diag(sig_y.0_sq[i], ncol(J))
        
        vcov.pst <- transform(.y = vcov.pst, trans = object@trans$scale)
        vcov.pst <- transform(.y = vcov.pst, trans = object@trans$scale)
        
        if (error %in% c("sd", "se"))
          return(list(fit = mu.pst, se = sqrt(diag(vcov.pst))))
        
        if (error %in% c("var", "cov", "vcov"))
          return(list(fit = mu.pst, vcov = vcov.pst))
        
        list(fit = mu.pst, se = sqrt(diag(vcov.pst)))
      })
    }
    
    if (length(ind.1) > 0) {
      KJ <- matrix(NA_real_, ncol(object@data$xt.1), ncol(xt.1))
      J  <- matrix(NA_real_, ncol(xt.1), ncol(xt.1))
      if (is.null(object@env$sampler.1)) sample.gpci.independent(object)
      pars <- extract(object@env$sampler.1, object@env$sampler.1@model_pars[seq_len(length(object@env$sampler.1@model_pars) - 1L)])
      
      origPars <- t(sapply(seq_len(nrow(pars$transformedPars)), function(i) devEnv$transformPars.1(pars$transformedPars[i,])))
      
      sig_f.1_sq <- exp(origPars[,1])
      origPars <- as.matrix(origPars[,-1])
      sig_y.1_sq <- exp(pars$log_sig_y_sq)
      mu <- pars$mu
      
      pst.1 <- lapply(seq_len(nrow(origPars)), function(i) {
        scales <- exp(-0.5 * origPars[i,])
        .Call(C_npci_grouped_updateCovMatrix, object@K.1, object@data$xt.1, object@data$xt.1, scales, sig_f.1_sq)
        
        .Call(C_npci_grouped_updateCovMatrix, KJ, object@data$xt.1, xt.1, scales, sig_f.1_sq)
        
        .Call(C_npci_grouped_updateLeftFactor, object@L.1, object@K.1)
        LKJ <- solve(object@L.1, KJ)
        mu.pst <- transform(.y = as.vector(mu[i] + crossprod(LKJ, solve(object@L.1, object@data$y.1 - mu[i]))),
                            trans = object@trans$inverse, simplify = TRUE)
        
        if (is.null(error)) return(mu.pst)
        
        .Call(C_npci_grouped_updateCovMatrix, J, xt.1, xt.1, scales, sig_f.1_sq)
        vcov.pst <- sig_y.1_sq[i] * (J - crossprod(LKJ))
        if (identical(error, "ppd")) vcov.pst <- vcov.pst + diag(sig_y.1_sq[i], ncol(J))
        
        vcov.pst <- transform(.y = vcov.pst, trans = object@trans$scale)
        vcov.pst <- transform(.y = vcov.pst, trans = object@trans$scale)
        
        if (error %in% c("sd", "se"))
          return(list(fit = mu.pst, se = sqrt(diag(vcov.pst))))
        
        if (error %in% c("var", "cov", "vcov"))
          return(list(fit = mu.pst, vcov = vcov.pst))
        
        list(fit = mu.pst, se = sqrt(diag(vcov.pst)))
      })
    }
    
    n.samp <- sum(sapply(object@env$sampler.0@stan_args, function(x) x$iter - x$warmup)) 
        
    if (is.null(error)) {
      mu.pst <- matrix(NA_real_, n.samp, nrow(x))
      if (length(ind.0) > 0) mu.pst[,ind.0] <- unlist(pst.0)
      if (length(ind.1) > 0) mu.pst[,ind.1] <- unlist(pst.1)
      
      return(res)
    }
    
    mu.pst <- matrix(NA_real_, n.samp, nrow(x))
    if (length(ind.0) > 0) mu.pst[,ind.0] <- t(sapply(pst.0, function(pst.i) pst.i$fit))
    if (length(ind.1) > 0) mu.pst[,ind.1] <- t(sapply(pst.1, function(pst.i) pst.i$fit))
    
    if (error %in% c("sd", "se")) {
      se.pst <- matrix(NA_real_, n.samp, nrow(x))
      if (length(ind.0) > 0) se.pst[,ind.0] <- t(sapply(pst.0, function(pst.i) pst.i$se))
      if (length(ind.1) > 0) se.pst[,ind.1] <- t(sapply(pst.1, function(pst.i) pst.i$se))
      
      return(list(fit = mu.pst, se = se.pst))
    }
    if (error %in% c("var", "cov", "vcov")) {
      browser()
      vcov.pst.0 <- sapply(pst.0, function(pst.i) pst.i$vcov)
      vcov.pst.1 <- sapply(pst.1, function(pst.i) pst.i$vcov)
    }
    return(pst)
  }
  
  if (NCOL(xt.0) > 0) {
    origPars.0 <- devEnv$transformPars.0(pars[1,])
  
    sig_f.0_sq <- exp(origPars.0[1L])
    scales.0 <- exp(-0.5 * origPars.0[-1L])
    
    dev.0 <- object@deviance.0(pars)
    mu.0.hat       <- attr(dev.0, "beta")
    sig_y.0_sq.hat <- attr(dev.0, "sig_y_sq")
    
    
    .Call(C_npci_grouped_updateCovMatrix, object@K.0, object@data$xt.0, object@data$xt.0, scales.0, sig_f.0_sq)
    
    KJ <- matrix(NA_real_, ncol(object@data$xt.0), ncol(xt.0))
    .Call(C_npci_grouped_updateCovMatrix, KJ, object@data$xt.0, xt.0, scales.0, sig_f.0_sq)
    
    .Call(C_npci_grouped_updateLeftFactor, object@L.0, object@K.0)
    LKJ <- solve(object@L.0, KJ)
    mu.pst.0 <- 
      transform(.y = as.vector(mu.0.hat + crossprod(LKJ, solve(object@L.0, object@data$y.0 - mu.0.hat))),
                     trans = object@trans$inverse, simplify = TRUE)
    
    if (!is.null(error)) {
      J  <- matrix(NA_real_, ncol(xt.0), ncol(xt.0))
      .Call(C_npci_grouped_updateCovMatrix, J, xt.0, xt.0, scales.0, sig_f.0_sq)
      vcov.pst.0 <- sig_y.0_sq.hat * (J - crossprod(LKJ))
      
      if (error %in% c("sd", "se")) {
        se.pst.0 <-
          transform(.y = as.vector(sqrt(diag(vcov.pst.0))), trans = object@trans$scale, simplify = TRUE)
      } else if (error %in% c("var", "cov", "vcov")) {
        vcov.pst.0 <- transform(.y = vcov.pst, trans = object@trans$scale)
        vcov.pst.0 <- transform(.y = vcov.pst, trans = object@trans$scale)
      }
    }
  }
  
  if (NCOL(xt.1) > 0) {
    origPars.1 <- devEnv$transformPars.1(pars[2,])
  
    sig_f.1_sq <- exp(origPars.1[1L])
    scales.1 <- exp(-0.5 * origPars.1[-1L])
    
    dev.1 <- object@deviance.1(pars)
    mu.1.hat       <- attr(dev.1, "beta")
    sig_y.1_sq.hat <- attr(dev.1, "sig_y_sq")
    
    
    .Call(C_npci_grouped_updateCovMatrix, object@K.1, object@data$xt.1, object@data$xt.1, scales.1, sig_f.1_sq)
    
    KJ <- matrix(NA_real_, ncol(object@data$xt.1), ncol(xt.1))
    .Call(C_npci_grouped_updateCovMatrix, KJ, object@data$xt.1, xt.1, scales.1, sig_f.1_sq)
    
    .Call(C_npci_grouped_updateLeftFactor, object@L.1, object@K.1)
    LKJ <- solve(object@L.1, KJ)
    mu.pst.1 <- 
      transform(.y = as.vector(mu.1.hat + crossprod(LKJ, solve(object@L.1, object@data$y.1 - mu.1.hat))),
                     trans = object@trans$inverse, simplify = TRUE)
    
    if (!is.null(error)) {
      J  <- matrix(NA_real_, ncol(xt.1), ncol(xt.1))
      .Call(C_npci_grouped_updateCovMatrix, J, xt.1, xt.1, scales.1, sig_f.1_sq)
      vcov.pst.1 <- sig_y.1_sq.hat * (J - crossprod(LKJ))
      
      if (error %in% c("sd", "se")) {
        se.pst.1 <-
          transform(.y = as.vector(sqrt(diag(vcov.pst.1))), trans = object@trans$scale, simplify = TRUE)
      } else if (error %in% c("var", "cov", "vcov")) {
        vcov.pst.1 <- transform(.y = vcov.pst, trans = object@trans$scale)
        vcov.pst.1 <- transform(.y = vcov.pst, trans = object@trans$scale)
      }
    }
  }
  
  mu.pst <- rep(NA_real_, nrow(x))
  if (length(ind.0) > 0) mu.pst[ind.0] <- mu.pst.0
  if (length(ind.1) > 0) mu.pst[ind.1] <- mu.pst.1
  
  if (is.null(error)) return(mu.pst)
  
  if (error %in% c("sd", "se")) {
    se.pst <- rep(NA_real_, nrow(x))
    if (length(ind.0) > 0) se.pst[ind.0] <- se.pst.0
    if (length(ind.1) > 0) se.pst[ind.1] <- se.pst.1
    return(list(fit = mu.pst, se = se.pst))
  }
  
  if (error %in% c("var", "cov", "vcov")) {
    vcov.pst <- matrix(0, nrow(x), nrow(x))
    if (length(ind.0) > 0) vcov.pst[ind.0, ind.0] <- vcov.pst.0
    if (length(ind.1) > 0) vcov.pst[ind.1, ind.1] <- vcov.pst.1
    return(list(fit = mu.pst, vcov = vcov.pst))
  }

  mu.pst
}

independentStanCode <- "
functions {
  real computeCovariance(int p, vector x, vector x_, vector scales) {
    vector[p] delta;
    
    delta <- x - x_;
    return exp(-dot_self(delta ./ scales));
  }
}
data {
  int<lower = 1> numObservations;
  int<lower = 1> numPredictors;
  matrix[numObservations, numPredictors] x;
  vector[numObservations] y;
  
  matrix[numPredictors + 1, numPredictors + 1] parL;
  vector[numPredictors + 1] parM;
}
transformed data {
  vector[numPredictors] xr[numObservations];
  for (i in 1:numObservations) {
    for (j in 1:numPredictors) {
      xr[i][j] <- x[i,j];
    }
  }
}
parameters {
  vector[numPredictors + 1] transformedPars;
  real log_sig_y_sq;
  real mu;
}
model {
  matrix[numObservations, numObservations] Sigma;
  vector[numObservations] mu_vec;
  
  {
    vector[numPredictors + 1] temp;
    vector[numPredictors] covPars;
    real sig_f_sq;
    real sig_y_sq;
    
    temp <- parL * transformedPars + parM;
    sig_f_sq <- exp(temp[1]);
    for (i in 1:numPredictors) covPars[i] <- exp(-0.5 * temp[i + 1]);
    sig_y_sq <- exp(log_sig_y_sq);
        
    for (i in 1:(numObservations - 1)) {
      for (j in (i + 1):numObservations) {
        Sigma[i,j] <- sig_f_sq * sig_y_sq * computeCovariance(numPredictors, xr[i], xr[j], covPars);
        Sigma[j,i] <- Sigma[i,j];
      }
      Sigma[i,i] <- sig_y_sq * (sig_f_sq + 1.0);
      mu_vec[i] <- mu;
    }
    Sigma[numObservations, numObservations] <- sig_y_sq * (sig_f_sq + 1.0);
    mu_vec[numObservations] <- mu;
  }
  
  for (i in 1:numPredictors) transformedPars[i] ~ cauchy(0, 2);   
  log_sig_y_sq ~ cauchy(0, 5);
  
  mu ~ cauchy(0, 10);
  
  y ~ multi_normal(mu_vec, Sigma);
}"

