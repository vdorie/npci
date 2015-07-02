ci.fit <- function(y, x, z, method, ...) {
  fit <- switch(method,
                naive1 = fitNaive1(y, x, z),
                naive2 = fitNaive2(y, x, z),
                bart   = fitBart(y, x, z, ...))
}



ci.estimate <- function(y, x, z, method, estimand, ...) {
  matchedCall <- match.call()
  
  fitCall <- stripCallArguments(retargetCall(matchedCall, quote(ci.fit)), "estimand")
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
    
    numSamples <- if (is.null(matchedCall$ndpost)) formals(predict.bartSampler)$ndpost else NULL
    
    pred <- predict(fit, x.test, z.test, ndpost = numSamples, keeptrainfits = TRUE, summarize = FALSE)
    
    return((pred$train[train.ind,] - pred$test) * sign)
  
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
    
    n <- NROW(x.test) %/% 2
    lin.comb <- 2 * z.test - 1
    
    pred <- predict(fit, x.test, z.test, "vcov")
    
    te <- pred$fit[z.test] - pred$fit[!z.test]
    se  <- sqrt(crossprod(lin.comb, pred$vcov) %*% lin.comb)[1] / n
    
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
  }
}
