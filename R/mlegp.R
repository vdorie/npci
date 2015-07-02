predict.gp <- function(object, x, error = NULL)
{
  if (missing(x)) x <- object$X
  
  meanFunc <- function(fit, x) { rep(if (identical(fit$constantMean, 1)) fit$mu[1] else 0, NROW(x)) }
  covFunc <- function(fit, x, y) {
    x <- as.matrix(x); y <- as.matrix(y)
    fit$sig2 * sapply(seq_len(NROW(y)), function(i)
      sapply(seq_len(NROW(x)), function(j) { delta <- x[j,] - y[i,]; exp(-sum(fit$beta * delta^2)) } )
    )
  }
  
  mu <- if (identical(object$constantMean, 1)) object$mu else rep(0, object$numObs)
  nu <- meanFunc(object, x)
  
  ## object$invVarMatrix = solve(K + diag(object$nugget, nrow(K)))
  
  ## K  <- covFunc(object, object$X, object$X) # not used
  J  <- covFunc(object, x, x)
  KJ <- covFunc(object, object$X, x)
  
  VKJ <- crossprod(KJ, object$invVarMatrix)
    
  mu.pst  <- nu + VKJ %*% (object$Z - mu)
  
  if (is.null(error)) return(mu.pst)
    
  vcov.pst <- J - VKJ %*% KJ
  
  if (error %in% c("sd", "se"))
    return(list(fit = mu.pst, se = sqrt(diag(vcov.pst))))
  
  if (error %in% c("var", "cov", "vcov"))
    return(list(fit = mu.pst, vcov = vcov.pst))
  
  warning("unrecognized error option")
  
  result
}
