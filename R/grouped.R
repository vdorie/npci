setClass("gpci.grouped", contains = "gpci",
         slots = list(deviance = "function",
                      K = "matrix", 
                      L = "matrix"),
         validity = function(object) {
           if (length(object@data$y) != NROW(object@K) || length(object@data$y) != NROW(object@L)) return("dimensions of cov matrix and factor must match length of 'y'")
           TRUE
         }
)


setMethod("initialize", "gpci.grouped",
  function(.Object, y, x, z)
{
  df <- as.data.frame(x)
  df$.y <- y
  trans <- getTransformations(df)
  df <- transform(df, trans$standardize)
  
  x <- as.matrix(df[,names(df) %not_in% ".y"])
  xt <- t(cbind(x, z = z))
  n <- length(y)
  
  .Object@data <- list(y = df$.y,
                       z = z,
                       x = x,
                       xt = xt,
                       x.mean = cbind(z, 1 - z))
  .Object@K <- matrix(NA_real_, n, n)
  .Object@L <- matrix(NA_real_, n, n)
  .Object@trans <- trans
  .Object@env <- new.env(parent = baseenv())
  
  devianceEnv <- args2env(baseenv(), object = .Object)
  devianceEnv$transformPars <- function(x) x
  
  .Object@deviance <- deviance.gpci.grouped
  environment(.Object@deviance) <- devianceEnv
  
  validObject(.Object)
  .Object
})

deviance.gpci.grouped <- function(pars)
{
  pars <- transformPars(pars)
  sig_f_sq <- exp(pars[1L])
  scales <- exp(-0.5 * pars[-1L])
  
  n <- length(object@data$y)
  .Call(npci:::C_npci_grouped_updateCovMatrix, object@K, object@data$xt, object@data$xt, scales, sig_f_sq)
  
  cholResult <- .Call(npci:::C_npci_grouped_updateLeftFactor, object@L, object@K)
  if (cholResult != 0L) return(.Machine$double.xmax * .Machine$double.eps^2)
  
  #(X'Sigma^-1X)^-1 X'Sigma^-1y
  #Sigma = LL', Sigma^-1 = L'^-1 L^-1
  #(X'L'^-1 L^-1 X)^-1 X'L'^-1 L^-1 y
  xr <- solve(object@L, object@data$x.mean)
  xtx <- crossprod(xr)
  
  beta.hat <- solve(xtx, crossprod(xr, solve(object@L, object@data$y)))
  mu <- as.vector(object@data$x.mean %*% beta.hat)
  # (y - x beta)' L'^-1 L^-1 (y - x beta)
  sig_y_sq.hat <- crossprod(solve(object@L, object@data$y - mu))[1] / n
  
  result <- n * (log(sig_y_sq.hat) + log(2 * pi)) + 2 * sum(log(diag(object@L))) + n
  attr(result, "beta")     <- beta.hat
  attr(result, "sig_y_sq") <- sig_y_sq.hat
  result
}

optimize.gpci.grouped <- function(object, n.testPoints = c(50L, 15L), n.modes = 4L, verbose = FALSE, start.only = FALSE, transform.pars = TRUE)
{
  if (is.null(object@env$opt)) object@env$opt <- list()
  if (length(object@env$opt) >= n.modes) return(invisible(NULL))
  
  n.covPars <- 1L + ncol(object@data$x)
  
  if (is.null(object@env$optStart)) {
    n.testPoints <- n.testPoints[1] + n.testPoints[2L] * n.covPars
    testMins  <- c(0, rep(-3, n.covPars))
    testMaxes <- c(10, rep(3, n.covPars))  
    
    testPoints <- matrix(runif(n.testPoints * (n.covPars + 1L), testMins, testMaxes), n.testPoints, byrow = TRUE)
    
    if (verbose) cat("calculating deviances\n")
    testDeviances <- t(sapply(seq_len(nrow(testPoints)), function(i) {
      object@deviance(testPoints[i,])
    }))
    
    devOrder <- order(testDeviances)
    testPoints    <- testPoints[devOrder,]
    testDeviances <- testDeviances[devOrder]
    
    n.testPoints <- ceiling(0.4 * n.testPoints)
    testPoints    <- testPoints[seq_len(n.testPoints),]
    testDeviances <- testDeviances[seq_len(n.testPoints)]
    
    kfit <- kmeans(testPoints, n.modes, nstart = 5L)
    object@env$optStart <- kfit$centers
  }
  
  if (identical(start.only, TRUE)) return(invisible(NULL))
  
  startTime <- Sys.time()
  for (i in seq.int(length(object@env$opt) + 1L, n.modes)) {
    if (verbose) cat("optimizing iteration ", i, "\n", sep = "")
    object@env$opt[[length(object@env$opt) + 1L]] <- optim(object@env$optStart[i,], object@deviance, method = "BFGS", hessian = TRUE)
  }
  endTime <- Sys.time()
  if (verbose) print(endTime - startTime)
  
  devs <- sapply(object@env$opt, function(opt) opt$value)
  object@env$opt <- object@env$opt[order(devs)]
  
  if (!identical(transform.pars, FALSE)) {
    devEnv <- environment(object@deviance)
    if (is.null(environment(devEnv$transformPars)$L)) {
      tryResult <- tryCatch(L <- solve(chol(object@env$opt[[1]]$hessian)), error = function(e) e)
      if (is(tryResult, "error")) L <- diag(1, n.covPars + 1L)
      
      x.0 <- object@env$opt[[1]]$par
      devEnv$transformPars <- function(x) L %*% x + x.0
      environment(devEnv$transformPars) <- args2env(baseenv(), L, x.0)
      for (i in seq_along(object@env$opt))
        object@env$opt[[i]]$par <- solve(L, object@env$opt[[i]]$par - x.0)
    }
  }
  
  invisible(NULL)
} 


createDevianceWrapper.gpci.grouped <- function(object, par.no, pars = object@env$opt[[1]]$par)
{
  wrapper <- function(x) {
    sapply(x, function(x.i) {
      pars.i <- pars
      pars.i[par.no] <- x.i
      deviance(pars.i)
    })
  }
  environment(wrapper) <- args2env(baseenv(), pars = pars, par.no = par.no, deviance = object@deviance)
  wrapper
}

sample.gpci.grouped <- function(object, n.chains = 4L, n.iter = 300L, n.burn = 150L, verbose = FALSE) {
  startingPoints <- "random"
  if (!is.null(object@env$opt) && length(object@env$opt) > 0) {
    startingPoints <- lapply(seq_along(object@env$opt), function(i) {
      temp <- object@deviance(object@env$opt[[i]]$par)
      list(transformedPars = object@env$opt[[i]]$par,
           log_sig_y_sq = log(attr(temp, "sig_y_sq")),
           mu_0 = attr(temp, "beta")[1L],
           mu_1 = attr(temp, "beta")[2L])
    })
    startingPoints <- rep_len(startingPoints, n.chains)
  }
  fit <- if (!is.null(object@env$sampler)) object@env$sampler else NA
  
  devEnv <- environment(object@deviance)
  
  .tempEnv <- NULL
  if (!exists("cpp_object_initializer", envir = .GlobalEnv)) {
    .tempEnv <- npci:::args2env(.GlobalEnv, cpp_object_initializer = Rcpp::cpp_object_initializer)
    attach(.tempEnv)
  }
  object@env$sampler <-
    stan(model_code = groupedStanCode,
         data = list(numObservations = length(object@data$y),
                     numPredictors   = ncol(object@data$x),
                     x = object@data$x,
                     y = object@data$y,
                     z = object@data$z,
                     parL = environment(devEnv$transformPars)$L,
                     parM = environment(devEnv$transformPars)$x.0),
         fit = fit,
         init = startingPoints,
         chains = n.chains,
         iter = n.iter,
         warmup = n.burn,
         verbose = verbose)
  
  if (!is.null(.tempEnv)) detach(.tempEnv)
  
  invisible(NULL)
}

predict.gpci.grouped <- function(object, x, z, error = NULL, pars = "sampler")
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
  xt <- t(cbind(x, z))
  x.mean <- cbind(z, 1 - z)
  
  devEnv <- environment(object@deviance)
  if (is.null(object@env$opt) && is.null(object@env$sampler)) optimize.gpci.grouped(object)
  
  if (identical(pars, "sampler")) {
    KJ <- matrix(NA_real_, ncol(object@data$xt), ncol(xt))
    J  <- matrix(NA_real_, ncol(xt), ncol(xt))
    if (is.null(object@env$sampler)) sample.gpci.grouped(object)
    pars <- extract(object@env$sampler, object@env$sampler@model_pars[seq_len(length(object@env$sampler@model_pars) - 1L)])
    
    origPars <- t(sapply(seq_len(nrow(pars$transformedPars)), function(i) devEnv$transformPars(pars$transformedPars[i,])))
    
    sig_f_sq <- exp(origPars[,1])
    origPars <- origPars[,-1]
    sig_y_sq <- exp(pars$log_sig_y_sq)
    mu_0 <- pars$mu_0
    mu_1 <- pars$mu_1
    
    pst <- lapply(seq_len(nrow(origPars)), function(i) {
      mu <- c(mu_1[i], mu_0[i])[object@data$z + 1]
      nu <- c(mu_1[i], mu_0[i])[z + 1]
      
      scales <- exp(-0.5 * origPars[i,])
      .Call(C_npci_grouped_updateCovMatrix, object@K, object@data$xt, object@data$xt, scales, sig_f_sq)
      
      .Call(C_npci_grouped_updateCovMatrix, KJ, object@data$xt, xt, scales, sig_f_sq)
      
      .Call(C_npci_grouped_updateLeftFactor, object@L, object@K)
      LKJ <- solve(object@L, KJ)
      mu.pst <- transform(.y = as.vector(nu + crossprod(LKJ, solve(object@L, object@data$y - mu))),
                          trans = object@trans$inverse, simplify = TRUE)
      
      if (is.null(error)) return(mu.pst)
      
      .Call(C_npci_grouped_updateCovMatrix, J, xt, xt, scales, sig_f_sq)
      vcov.pst <- sig_y_sq[i] * (J - crossprod(LKJ))
      if (identical(error, "ppd")) vcov.pst <- vcov.pst + diag(sig_y_sq[i], ncol(J))
      
      vcov.pst <- transform(.y = vcov.pst, trans = object@trans$scale)
      vcov.pst <- transform(.y = vcov.pst, trans = object@trans$scale)
      
      if (error %in% c("sd", "se"))
        return(list(fit = mu.pst, se = sqrt(diag(vcov.pst))))
      
      if (error %in% c("var", "cov", "vcov"))
        return(list(fit = mu.pst, vcov = vcov.pst))
      
      list(fit = mu.pst, se = sqrt(diag(vcov.pst)))
    })
    
    if (is.null(error)) return(unlist(pst))
    
    mu.pst <- t(sapply(pst, function(pst.i) pst.i$fit))
    if (error %in% c("sd", "se")) {
      se.pst <- t(sapply(pst, function(pst.i) pst.i$se))
      return(list(fit = mu.pst, se = se.pst))
    }
    if (error %in% c("var", "cov", "vcov")) {
      browser()
      vcov.pst <- sapply(pst, function(pst.i) pst.i$vcov)
    }
    return(pst)
  }
  
  origPars <- devEnv$transformPars(pars)
      

  sig_f_sq <- exp(origPars[1L]) ## deviance may use a whitening transformation
  scales <- exp(-0.5 * origPars[-1L])

  dev <- object@deviance(pars)
  beta.hat     <- attr(dev, "beta")
  sig_y_sq.hat <- attr(dev, "sig_y_sq")
  
  mu <- object@data$x.mean %*% beta.hat
  nu <- x.mean %*% beta.hat
  
  .Call(C_npci_grouped_updateCovMatrix, object@K, object@data$xt, object@data$xt, scales, sig_f_sq)
  
  KJ <- matrix(NA_real_, ncol(object@data$xt), ncol(xt))
  .Call(C_npci_grouped_updateCovMatrix, KJ, object@data$xt, xt, scales, sig_f_sq)
  
  
  .Call(C_npci_grouped_updateLeftFactor, object@L, object@K)
  LKJ <- solve(object@L, KJ)
  mu.pst <- 
    transform(.y = as.vector(nu + crossprod(LKJ, solve(object@L, object@data$y - mu))),
              trans = object@trans$inverse, simplify = TRUE)
  
  if (is.null(error)) return(mu.pst)
  
  J  <- matrix(NA_real_, ncol(xt), ncol(xt))
  .Call(C_npci_grouped_updateCovMatrix, J, xt, xt, scales, sig_f_sq)
  vcov.pst <- sig_y_sq.hat * (J - crossprod(LKJ))
  
  if (error %in% c("sd", "se")) {
    se.pst <-
      transform(.y = as.vector(sqrt(diag(vcov.pst))), trans = object@trans$scale, simplify = TRUE)
    return(list(fit = mu.pst, se = se.pst))
  }
  
  if (error %in% c("var", "cov", "vcov")) {
    vcov.pst <-
      transform(.y = vcov.pst, trans = object@trans$scale)
    vcov.pst <-
      transform(.y = vcov.pst, trans = object@trans$scale)
    return(list(fit = mu.pst, vcov = vcov.pst))
  }
  
  mu.pst
}

groupedStanCode <- "
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
  vector[numObservations] z;
  
  matrix[numPredictors + 2, numPredictors + 2] parL;
  vector[numPredictors + 2] parM;
}
transformed data {
  int<lower = 1> numPredictors_gp;
  vector[numPredictors + 1] xr[numObservations]; // x by row, with z added as a column
  
  numPredictors_gp <- numPredictors + 1;
  for (i in 1:numObservations) {
    for (j in 1:numPredictors) {
      xr[i][j] <- x[i,j];
    }
    xr[i][numPredictors_gp] <- z[i];
  }
}
parameters {
  vector[numPredictors_gp + 1] transformedPars;
  real log_sig_y_sq;
  real mu_0;
  real mu_1;
}
model {
  matrix[numObservations, numObservations] Sigma;
  vector[numObservations] mu;
  
  {
    vector[numPredictors_gp + 1] temp;
    vector[numPredictors_gp] covPars;
    real sig_f_sq;
    real sig_y_sq;
    
    temp <- parL * transformedPars + parM;
    sig_f_sq <- exp(temp[1]);
    for (i in 1:numPredictors_gp) covPars[i] <- exp(-0.5 * temp[i + 1]);
    sig_y_sq <- exp(log_sig_y_sq);
    // print(\"sig_f_sq: \", sig_f_sq, \", sig_y_sq: \", sig_y_sq, \", covPars: \", covPars);
    // print(\"x_1: \", xr[1]);
    
    for (i in 1:(numObservations - 1)) {
      for (j in (i + 1):numObservations) {
        Sigma[i,j] <- sig_f_sq * sig_y_sq * computeCovariance(numPredictors_gp, xr[i], xr[j], covPars);
        Sigma[j,i] <- Sigma[i,j];
      }
      Sigma[i,i] <- sig_y_sq * (sig_f_sq + 1.0);
      mu[i] <- if_else(z[i] != 0.0, mu_1, mu_0);
    }
    Sigma[numObservations, numObservations] <- sig_y_sq * (sig_f_sq + 1.0);
    mu[numObservations] <- if_else(z[numObservations] != 0.0, mu_1, mu_0);
    
    // print(\"sigma: \", Sigma);
  }
  
  for (i in 1:numPredictors_gp) transformedPars[i] ~ cauchy(0, 2);   
  log_sig_y_sq ~ cauchy(0, 5);
  
  mu_0 ~ cauchy(0, 10);
  mu_1 ~ cauchy(0, 10);
  
  y ~ multi_normal(mu, Sigma);
}"

