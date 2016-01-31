fitBTGP <- function(y, x, z, nskip = 200L, ndpost = 200L, nthin = 5L) {
  df <- as.data.frame(x)
  groupMeans <- c(mean(y[z == 0]), mean(y[z == 1]))
  df$.y <- y - groupMeans[z + 1]
  trans <- getTransformations(df)
  df <- transform(df, trans$standardize)
  
  y <- df$.y
  x <- cbind(as.matrix(df[,names(df) %not_in% ".y"]), z)
  BTE <- c(nskip * nthin, (nskip + ndpost) * nthin, nthin)
  
  result <- namedList(x, y, groupMeans, trans, BTE)
  class(result) <- "btgp"
  
  result
}

predict.btgp <- function(object, x, z, pred.n = TRUE)
{
  x <- as.matrix(x)
  if (is.null(colnames(x)))
    colnames(x) <- names(object$trans$standardize)[names(object$trans$standardize) %not_in% ".y"]
  
  x <- as.matrix(transform(as.data.frame(x), object$trans$standardize))
  z <- rep_len(z, nrow(x))
  x <- cbind(x, z)
  
  fit <- btgp(object$x, object$y, x, BTE = object$BTE, zcov = TRUE, pred.n = pred.n, krige = FALSE,
              verb = 0L, m0r1 = FALSE, meanfn = "constant")
  
  mu <- if (pred.n) c(fit$Zp.mean, fit$ZZ.mean) else fit$ZZ.mean
  mu <- transform(.y = mu, trans = object$trans$inverse, simplify = TRUE)
  mu <- mu + object$groupMeans[(if (pred.n) c(object$x[,ncol(object$x)], z) else z) + 1]
  
  n <- nrow(object$x)
  m <- nrow(x)
  extent <- m + if (pred.n) n else 0
  Sig <- matrix(NA, extent, extent)
  if (pred.n) {
    Sig[seq_len(n), seq_len(n)] <- fit$Zp.s2
    Sig[seq_len(n), n + seq_len(m)] <- fit$ZpZZ.s2
    Sig[n + seq_len(m), seq_len(n)] <- t(fit$ZpZZ.s2)
    Sig[n + seq_len(m), n + seq_len(m)] <- fit$ZZ.s2
  } else {
    Sig[seq_len(m), seq_len(m)] <- fit$ZZ.s2
  }
  Sig <-
      transform(.y = Sig, trans = object$trans$scale)
    vcov.pst <-
      transform(.y = Sig, trans = object$trans$scale)
  
  list(fit = mu, vcov = Sig)
}
