"%not_in%" <- function(x, table) match(x, table, nomatch = 0) == 0

namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

retargetCall <- function(call, symbol) {
  call[[1]] <- symbol
  call
}

stripCallArguments <- function(call, ...) {
  if (missing(call)) stop("call cannot be missing")
  extraArguments <- as.character(list(...))
  if (length(extraArguments) == 0) return(call)
  
  call[names(call) %not_in% extraArguments]
}

substituteTermInLanguage <- function(lang, oldTerm, newTerm)
{
  for (i in seq_along(lang)) {
    if (is.symbol(lang[[i]])) {
      if (lang[[i]] == oldTerm) lang[[i]] <- newTerm
    } else if (is.language(lang[[i]])) {
      lang[[i]] <- substituteTermInLanguage(lang[[i]], oldTerm, newTerm)
    }
  }
  lang
}

args2env <- function(parent, ...)
{
  parList <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(parList))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  names(parList) <- resultNames

  list2env(parList, parent = parent)
}
