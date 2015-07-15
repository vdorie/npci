fitNaive2 <- function(y, x, z) {
  df <- as.data.frame(x)
  groupMeans <- c(mean(y[z == 0]), mean(y[z == 1]))
  df$y <- y - groupMeans[z + 1]
  
  df.0 <- subset(df, z == 0)
  df.1 <- subset(df, z == 1)
  
  trans.0 <- getTransformations(df.0)
  trans.1 <- getTransformations(df.1)
  
  df.0 <- transform(df.0, trans.0$standardize)
  df.1 <- transform(df.1, trans.1$standardize)
  
  sig.est.0 <- summary(lm(y ~ ., df.0))$sigma
  sig.est.1 <- summary(lm(y ~ ., df.1))$sigma
  
  covCols <- names(df) %not_in% "y"
  
  ignored <- NULL
  stringConnection <- textConnection("ignored", "w", local = TRUE)
  sink(stringConnection)
  
  fit.0 <- mlegp(df.0[,covCols], df.0$y, nugget = sig.est.0^2, verbose = 0)
  fit.1 <- mlegp(df.1[,covCols], df.1$y, nugget = sig.est.1^2, verbose = 0)
  
  sink()
  close(stringConnection)
  
  result <- namedList(fit.0, fit.1, trans.0, trans.1, groupMeans, z)
  class(result) <- "naive2"
  result
}

predict.naive2 <- function(object, x, z, error = NULL, ...) {
  fit.0   <- object$fit.0
  fit.1   <- object$fit.1
  trans.0 <- object$trans.0
  trans.1 <- object$trans.1
  groupMeans <- object$groupMeans
  
  df <- as.data.frame(x)
  n <- nrow(df)
  z <- rep_len(z, n)
  
  zeroRows <- which(z == 0)
  oneRows  <- which(z == 1)
  
  result <- numeric(n)
  if (!is.null(error) && error %in% c("sd", "se"))
    { error <- "se"; se <- numeric(n) }
  if (!is.null(error) && error %in% c("var", "cov", "vcov"))
    { error <- "vcov"; vcov <- matrix(0, n, n) }
  
  if (length(zeroRows) > 0) {
    df.0 <- transform(subset(df, z == 0), trans.0$standardize)
    pred <- predict(fit.0, df.0, error = error)
    
    y.hat <- if (!is.list(pred)) pred else pred$fit
    
    result[zeroRows] <- transform(y = y.hat, trans = trans.0$inverse, simplify = TRUE) + groupMeans[1]
    
    if (identical(error, "se"))
      se[zeroRows] <- transform(y = pred$se, trans = trans.0$scale, simplify = TRUE)
    else if (identical(error, "vcov"))
      vcov[zeroRows, zeroRows] <- transform(y = transform(y = pred$vcov, trans = trans.0$scale), trans = trans.0$scale)
  }
  if (length(oneRows)  > 0) {
    df.1 <- transform(subset(df, z == 1),  trans.1$standardize)
    pred <- predict(fit.1, df.1, error = error)
    
    y.hat <- if (!is.list(pred)) pred else pred$fit
    
    result[oneRows] <- transform(y = y.hat, trans = trans.1$inverse, simplify = TRUE) + groupMeans[2]
    
    if (identical(error, "se"))
      se[oneRows] <- transform(y = pred$se, trans = trans.1$scale, simplify = TRUE)
    else if (identical(error, "vcov"))
      vcov[oneRows, oneRows] <- transform(y = transform(y = pred$vcov, trans = trans.1$scale), trans = trans.1$scale)
  }
  
  if (is.null(error)) return(result)
  if (identical(error, "se")) return(namedList(fit = result, se))
  if (identical(error, "vcov")) return(namedList(fit = result, vcov))
  
  result
}

fitted.naive2 <- function(object, se = FALSE, ...) {
  fit.0   <- object$fit.0
  fit.1   <- object$fit.1
  trans.0 <- object$trans.0
  trans.1 <- object$trans.1
  groupMeans <- object$groupMeans
  z <- object$z
  
  n <- length(z)
    
  zeroRows <- which(z == 0)
  oneRows  <- which(z == 1)

  result <- if (identical(se, TRUE)) list(fit = numeric(n), se.fit = numeric(n)) else numeric(n)
  
  if (length(zeroRows) > 0) {
    temp <- predict(fit.0, se.fit = se)
  
    if (identical(se, TRUE)) {
      result$fit[zeroRows]    <- transform(y = temp$fit,    trans = trans.0$inverse, simplify = TRUE) + groupMeans[1]
      result$se.fit[zeroRows] <- transform(y = temp$se.fit, trans = trans.0$scale,   simplify = TRUE)
    } else {
      result[zeroRows] <- transform(y = temp, trans = trans.0$inverse, simplify = TRUE) + groupMeans[1]
    }
  }
  if (length(oneRows)  > 0) {
    temp <- predict(fit.1, se.fit = se)
    if (identical(se, TRUE)) {
      result$fit[oneRows]     <- transform(y = temp$fit,    trans = trans.1$inverse, simplify = TRUE) + groupMeans[2]
      result$se.fit[oneRows]  <- transform(y = temp$se.fit, trans = trans.1$scale,   simplify = TRUE)
    } else {
      result[oneRows]  <- transform(y = temp, trans = trans.1$inverse, simplify = TRUE) + groupMeans[2]
    }
  }
  result
}
