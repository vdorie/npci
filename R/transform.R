getTransformations <- function(x, ignore = NULL, type = "norm")
{
  ## cc - continuous columns
  cc <- sapply(x, function(col) is.numeric(col) && length(unique(col)) > 2)
  if (!is.null(ignore)) cc <- cc & names(x) %not_in% ignore
  
  if (identical(type, "norm")) {
    fnEnvs <- lapply(names(x)[cc], function(name) { args2env(baseenv(), mu = mean(x[[name]]), sig = sd(x[[name]])) })
    names(fnEnvs) <- names(x)[cc]
    standardize <- lapply(names(x)[cc], function(name) {
      res <- function(x) { (x - mu) / sig }
      environment(res) <- fnEnvs[[name]]
      res
    })
    inverse <- lapply(names(x)[cc], function(name) {
      res <- function(y) { y * sig + mu }
      environment(res) <- fnEnvs[[name]]
      res
    })
    scale <- lapply(names(x)[cc], function(name) {
      res <- function(y) y * sig
      environment(res) <- fnEnvs[[name]]
      res
    })
    names(standardize) <- names(inverse) <- names(scale) <- names(x)[cc]
  } else if (identical(type, "cube")) {
    fnEnvs <- lapply(seq_along(x), function(j) { args2env(baseenv(), m = min(x[[j]]), Mm = max(x[[j]]) - min(x[[j]])) })
    standardize <- lapply(seq_along(x), function(j) {
      res <- function(x) { (x - m) / Mm }
      environment(res) <- fnEnvs[[j]]
      res
    })
    inverse <- lapply(seq_along(x), function(j) {
      res <- function(y) { y * Mm + m }
      environment(res) <- fnEnvs[[j]]
      res
    })
    scale <- lapply(seq_along(x), function(j) {
      res <- function(y) { y * Mm }
      environment(res) <- fnEnvs[[j]]
      res
    })
    names(standardize) <- names(inverse) <- names(scale) <- names(x)
  } else {
    stop("unrecognized transformation type: ", type)
  }

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
