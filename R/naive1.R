fitNaive1 <- function(y, x, z) {
  df <- as.data.frame(x)
  groupMeans <- c(mean(y[z == 0]), mean(y[z == 1]))
  df$y <- y - groupMeans[z + 1]
  
  trans <- getTransformations(df)
  df <- transform(df, trans$standardize)
  
  df$z <- z
  
  sig.est <- summary(lm(y ~ ., df))$sigma
  
  covCols <- names(df) %not_in% "y"
  
  ignored <- NULL
  stringConnection <- textConnection("ignored", "w", local = TRUE)
  sink(stringConnection)
  
  fit <- mlegp(df[, covCols], df$y, nugget = sig.est^2, verbose = 0)
  
  sink()
  close(stringConnection)
  
  result <- namedList(fit, trans, groupMeans)
  class(result) <- "naive1"
  result
}

predict.naive1 <- function(object, x, z, error = NULL, ...) {
  fit   <- object$fit
  trans <- object$trans
  groupMeans <- object$groupMeans
  
  df <- as.data.frame(x)
  df <- transform(df, trans$standardize)
  df$z <- rep_len(z, nrow(df))
  
  x.test <- as.matrix(df)
  
  pred <- predict(fit, x.test, error)
  
  y.hat <- if (!is.list(pred)) pred else pred$fit
  y.hat <- transform(y = y.hat, trans = trans$inverse, simplify = TRUE) + groupMeans[z + 1]
  
  if (!is.list(pred)) return(y.hat)
  
  pred$fit <- y.hat
  
  if (names(pred)[[2]] == "se") {
    pred$se <- transform(y = pred$se, trans = trans$scale, simplify = TRUE)
  } else if (names(pred)[[2]] == "vcov") {
    pred$vcov <- transform(y = pred$vcov, trans = trans$scale)
    pred$vcov <- transform(y = pred$vcov, trans = trans$scale)
  }
  
  pred
}

fitted.naive1 <- function(object, se = FALSE, ...) {
  fit <- object$fit
  trans <- object$trans
  groupMeans <- object$groupMeans
  
  y.hat <- predict(fit, se.fit = se)
  offset <- groupMeans[fit$X[, "z"]]
  
  if (identical(se, TRUE))
    return(list(fit    = transform(y = y.hat$fit,    trans = trans$inverse, simplify = TRUE) + offset,
                se.fit = transform(y = y.hat$se.fit, trans = trans$scale,   simplify = TRUE)))
  
  transform(y = y.hat, trans = trans$inverse, simplify = TRUE) + offset
}
