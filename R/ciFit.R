ci.fit <- function(y, x, z, method, estimand, prob.z = NULL, ...) {
  if (missing(estimand)) estimand <- NULL
  else if (!is.null(estimand) && estimand %not_in% c("ate", "att", "atc"))
    stop("estimand '", estimand, "' not recognized")
  
  w.est <- if (is.null(prob.z) || is.null(estimand)) NULL else
    switch(estimand,
           ate = ifelse(z, 1 / prob.z, 1 / (1 - prob.z)),
           att = ifelse(z, 1, prob.z / (1 - prob.z)),
           atc = ifelse(z, (1 - prob.z) / prob.z, 1))

  fit <- switch(method,
                naive1 = fitNaive1(y, x, z),
                naive2 = fitNaive2(y, x, z),
                bart   = fitBart(y, x, z, weights = w.est, ...),
                grouped = {
                  verbose <- if (is.null(verbose <- list(...)$verbose)) FALSE else verbose
                  m <- new("gpci.grouped", y, x, z)
                  optimize.gpci.grouped(m, verbose = verbose)
                  # sample.gpci.grouped(m, verbose = verbose)
                  m
               },
               independent = {
                 verbose <- if (is.null(verbose <- list(...)$verbose)) FALSE else verbose
                 m <- new("gpci.independent", y, x, z)
                 optimize.gpci.independent(m, verbose = verbose)
                  # sample.gpci.independent(m, verbose = verbose)
                 m
               },
               btgp = fitBTGP(y, x, z, ...)
  )
}



ci.estimate <- function(y, x, z, method, estimand, prob.z = NULL, ...) {
  matchedCall <- match.call()
  
  if (estimand %not_in% c("ate", "att", "atc"))
    stop("estimand '", estimand, "' not recognized")
  
  # fitCall <- stripCallArguments(retargetCall(matchedCall, quote(ci.fit)), "estimand")
  fitCall <- retargetCall(matchedCall, quote(ci.fit))
  callingEnv <- parent.frame(1L)
  fit <- eval(fitCall, callingEnv)
    
  if (identical(method, "bart")) {
    if (identical(estimand, "ate")) {
      train.ind <- seq_len(NROW(x))
      x.test <- x
      z.test <- 1 - z
      sign <- ifelse(z == 1, 1, -1) ## column vector gets reused when multiplied
    } else if (identical(estimand, "att")) {
      train.ind <- which(z == 1)
      x.test <- subset(x, z == 1)
      z.test <- 0
      sign <- 1
    } else if (identical(estimand, "atc")) {
      train.ind <- which(z == 0)
      x.test <- subset(x, z == 0)
      z.test <- 1
      sign <- -1
    }
    
    numSamples <- if (is.null(matchedCall$ndpost)) formals(predict.bartSampler)$ndpost else matchedCall$ndpost 
    
    pred <- predict(fit, x.test, z.test, ndpost = numSamples, keeptrainfits = TRUE, summarize = FALSE)
    return((pred$train[train.ind,] - pred$test) * sign)
  
  } else if (identical(method, "btgp")) {
    if (identical(estimand, "ate")) {
      train.ind <- seq_len(NROW(x))
      x.test <- x
      z.test <- 1 - z
      sign <- ifelse(z == 1, 1, -1)
      lin.comb <- 2 * z - 1
    } else if (identical(estimand, "att")) {
      train.ind <- which(z == 1)
      x.test <- subset(x, z == 1)
      z.test <- 0
      sign <- 1
      lin.comb <- z
    } else if (identical(estimand, "atc")) {
      train.ind <- which(z == 0)
      x.test <- subset(x, z == 0)
      z.test <- 1
      sign <- -1
      lin.comb <- 1 - z
    }
    test.ind <- NROW(x) + seq_len(NROW(x.test))
    lin.comb <- c(lin.comb, 2 * rep_len(z.test, NROW(x.test)) - 1)
    
    pred <- predict(fit, x.test, z.test)
    
    te <- (pred$fit[train.ind] - pred$fit[test.ind]) * sign
    se <- sqrt(crossprod(lin.comb, pred$vcov) %*% lin.comb)[1] / NROW(x.test)
    
    return(namedList(te, se))
    
  } else if (identical(method, "naive1")) {
    if (identical(estimand, "ate")) {
      x.test <- rbind(x, x)
      z.test <- c(z, 1 - z)
    } else if (identical(estimand, "att")) {
      temp <- subset(x, z == 1)
      x.test <- rbind(temp, temp)
      z.test <- c(rep(1, NROW(temp)), rep(0, NROW(temp)))
    } else if (identical(estimand, "atc")) {
      temp <- subset(x, z == 0)
      x.test <- rbind(temp, temp)
      z.test <- c(rep(1, NROW(temp)), rep(0, NROW(temp)))
    }
    
    n <- NROW(x.test) %/% 2L
    lin.comb <- 2 * z.test - 1
    
    pred <- predict(fit, x.test, z.test, "vcov")
    
    te <- pred$fit[z.test == 1] - pred$fit[z.test == 0]
    se <- sqrt(crossprod(lin.comb, pred$vcov) %*% lin.comb)[1] / n
    
    return(namedList(te, se))
    
  } else if (identical(method, "naive2")) {
    if (identical(estimand, "ate")) {
      x.test <- x
    } else if (identical(estimand, "att")) {
      x.test <- subset(x, z == 1)
    } else if (identical(estimand, "atc")) {
      x.test <- subset(x, z == 0)
    }
    
    n <- NROW(x.test)
    lin.comb <- rep(1, n)
    
    pred.0 <- predict(fit, x.test, 0, "vcov")
    pred.1 <- predict(fit, x.test, 1, "vcov")
    
    te <- pred.1$fit - pred.0$fit
    se <- sqrt(crossprod(lin.comb, pred.1$vcov) %*% lin.comb + crossprod(lin.comb, pred.0$vcov) %*% lin.comb)[1] / n
    
    return(namedList(te, se))
  } else if (identical(method, "grouped")) {
    if (identical(estimand, "ate")) {
      x.test <- rbind(x, x)
      z.test <- c(z, 1 - z)
    } else if (identical(estimand, "att")) {
      temp <- subset(x, z == 1)
      x.test <- rbind(temp, temp)
      z.test <- c(rep(1, NROW(temp)), rep(0, NROW(temp)))
    } else if (identical(estimand, "atc")) {
      temp <- subset(x, z == 0)
      x.test <- rbind(temp, temp)
      z.test <- c(rep(1, NROW(temp)), rep(0, NROW(temp)))
    }
    
    n <- NROW(x.test) %/% 2L
    lin.comb <- 2 * z.test - 1
    
    pred <- predict(fit, x.test, z.test, "vcov", pars = fit@env$opt[[1]]$par)
    
    te <- pred$fit[z.test == 1] - pred$fit[z.test == 0]
    v <- (crossprod(lin.comb, pred$vcov) %*% lin.comb)[1]
    if (v < 0.0) browser()
    se <- sqrt(v) / n
    
    return(namedList(te, se, fit))
  } else if (identical(method, "independent")) {
    if (identical(estimand, "ate")) {
      x.test <- x
    } else if (identical(estimand, "att")) {
      x.test <- subset(x, z == 1)
    } else if (identical(estimand, "atc")) {
      x.test <- subset(x, z == 0)
    }
    
    n <- NROW(x.test)
    lin.comb <- rep(1, n)
    
    pred.0 <- predict(fit, x.test, 0, "vcov", pars = rbind(fit@env$opt.0[[1]]$par, fit@env$opt.1[[1]]$par))
    pred.1 <- predict(fit, x.test, 1, "vcov", pars = rbind(fit@env$opt.0[[1]]$par, fit@env$opt.1[[1]]$par))
    
    te <- pred.1$fit - pred.0$fit
    se <- sqrt(crossprod(lin.comb, pred.1$vcov) %*% lin.comb + crossprod(lin.comb, pred.0$vcov) %*% lin.comb)[1] / n
    
    return(namedList(te, se, fit))
  }
}
