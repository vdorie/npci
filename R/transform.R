getTransformations <- function(x, ignore = NULL)
{
  ## cc - continuous columns
  cc <- sapply(x, function(col) is.numeric(col) && length(unique(col)) > 2)
  if (!is.null(ignore)) cc <- cc & names(x) %not_in% ignore
  
  standardize <- lapply(names(x)[cc], function(name) {
    res <- function(x) { (x - mu) / sig }
    environment(res) <- list2env(list(mu = mean(x[[name]]), sig = sd(x[[name]])), parent = baseenv())
    res
  })
  names(standardize) <- names(x)[cc]
  inverse <- lapply(names(x)[cc], function(name) {
    res <- function(y) { y * sig + mu }
    environment(res) <- list2env(list(mu = mean(x[[name]]), sig = sd(x[[name]])), parent = baseenv())
    res
  })
  names(inverse) <- names(x)[cc]
  scale <- lapply(names(x)[cc], function(name) {
    res <- function(x) x * sig
    environment(res) <- list2env(list(sig = sd(x[[name]])), parent = baseenv())
    res
  })
  names(scale) <- names(x)[cc]

  namedList(standardize, inverse, scale)
}

transform <- function(data, trans, simplify = FALSE, ...)
{
  if (missing(data)) {
    data <- list(...)
    if (length(data) == 1 && is.matrix(data[[1]]))
      return(trans[[names(data)[[1]]]](data[[1]]))
    
    lengths <- sapply(data, length)
    if (all(lengths == lengths[1])) data <- as.data.frame(data)
  }
    
  for (name in names(data))
    if (!is.null(trans[[name]])) data[[name]] <- trans[[name]](data[[name]])
  if (simplify) {
    if (NCOL(data) == 1) {
      if (is.list(data)) as.numeric(data[[1]]) else as.numeric(data)
    } else as.matrix(data)
  } else data
}
