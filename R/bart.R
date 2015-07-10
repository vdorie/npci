fitBart <- function(y, x, z, ntree = 50L, nskip = 200L) {
  df <- as.data.frame(x)
  groupMeans <- c(mean(y[z == 0]), mean(y[z == 1]))
  df$y <- y - groupMeans[z + 1]
  df$z <- z
  sampler <- dbarts(y ~ ., df, control = dbartsControl(n.trees = ntree, n.burn = nskip, updateState = FALSE))
  sampler$run(nskip, 0)
  
  result <- namedList(sampler, groupMeans)
  class(result) <- "bartSampler"
  result
}

predict.bartSampler <- function(object, x, z, ndpost = 400L, keeptrainfits = TRUE, 
                                summarize = TRUE)
{
  sampler <- object$sampler
  groupMeans <- object$groupMeans
  
  df <- as.data.frame(x)
  if (is.null(colnames(x))) colnames(df) <- colnames(sampler$data@x)[-ncol(sampler$data@x)]
  df$z <- rep_len(z, nrow(df))
  
  sampler$setTestPredictor(df)
  
  samples <- sampler$run(0, ndpost)
  
  testSamples <- samples$test + matrix(groupMeans[df$z + 1], nrow(df), ndpost)
  trainingSamples <-
    if (!identical(keeptrainfits, TRUE))
      NULL
    else
      samples$train + matrix(groupMeans[sampler$data@x[,"z"] + 1], nrow(sampler$data@x), ndpost)
  
  if (identical(summarize, FALSE) || is.null(summarize))
    return (if (is.null(trainingSamples)) testSamples else list(train = trainingSamples, test = testSamples))
  
  if (identical(summarize, TRUE)) {
     summarize <- list(mean  = quote(mean(x)),
                       bound = quote(quantile(x, c(0.025, 0.975))))
  } else if (is.numeric(summarize)) {
    s.call <- quote(quantile(x, q))
    s.call[[3]] <- summarize
    summarize <- list(mean  = quote(mean(x)),
                      bound = s.call)
  } else if (is.symbol(summarize)) {
    s.call <- quote(dummy(x))
    s.call[[1]] <- summarize
    summarize <- s.call
  } else if (is.function(summarize)) {
    s.call <- quote(dummy(x))
    s.call[[1]] <- match.call()$summarize
    summarize <- s.call
  }
    
  evaluator <- function(x) {
    f.x <- substituteTermInLanguage(f, quote(x), x)
    eval(f.x)
  }
  environment(evaluator) <- new.env()
  functionEnv <- environment(evaluator)
  
  if (is.list(summarize) && length(summarize) == 1) summarize <- summarize[[1]]
  
  if (!is.list(summarize)) {
    functionEnv$f <- summarize
    result <- apply(testSamples, 1, evaluator)
    if (!is.null(trainingSamples))
      result <- list(test  = result,
                     train = apply(trainingSamples, 1, evaluator))
    return(result)
  }
  
  result <- vector("list", length(summarize))
  names(result) <- names(summarize)
  for (i in seq_along(summarize)) {
    functionEnv$f <- summarize[[i]]

    result.i <- apply(testSamples, 1, evaluator)
    if (!is.null(trainingSamples)) {
      result.i <- list(test  = result,
                       train = apply(trainingSamples, 1, evaluator))
    }
    
    result[[i]] <- result.i
  }
  result
}
